% Make figure 3 of the attribution manuscript, showing: 
% (a) the bed elevation
% (b) the ice geometry in stages 1 and 2 of the calibration
% (c) forcing profiles
% (d) zoom in of the ice shelf geometry from (b)
% (e) mean melt rate as a function of M

%
% Alex Bradley, 10/11/22. MIT License.

%
% Preliminaries
%
addpath('plottools')

%
% Constants
%
rhoi = 918;  %ice density
rhow = 1028; %water density

%% Fetch the data
%
% Get the results for the no-melt stage of the calibration
%
fpath = strcat('data/ATTR_00000/outfile.nc');
bed   = ncread(fpath, 'b', [1, 1, 1], [Inf, Inf, 1]); %bed topo
h     = ncread(fpath, 'h', [1, 1, 1], [Inf,Inf,Inf]); %ice thickness
h     = squeeze(h(:,:,end)); %ice thickness at the final timestep (in steady state)
gr    = ncread(fpath, 'grfrac', [1,1,1],[Inf,Inf,Inf]); %grounded fraction at all times in simulation
gr    = squeeze(gr(:,:,end)); %grounded fraction at the final timestep (in steady state)
s     = ncread(fpath, 's', [1,1,1],[Inf,Inf,Inf]); %surface elevation through the simulation
s     = squeeze(s(:,:,end)); %ice surface at the final timestep (in steady stat
b     = s - h; %ice bottom
x     = ncread(fpath, 'x'); %x co-ordinates
y     = ncread(fpath, 'y'); %x co-ordinates

%
% Get the results for post calibration runs
%
run_nums = ["08002", "09002", "10002", "11002", "12002"];
sz       = size(run_nums);
hf       = cell(sz);  %ice thickness
sf       = cell(sz);  %ice surface
mf       = cell(sz);  %melt rate
grf      = cell(sz);  %grounded fraction
bf 	     = cell(sz);  %ice bottom
for i = 1:sz(2)
	%filename
	fname =  strcat('data/ATTR_', run_nums(i), '/outfile.nc');

	%get thickness, surface, grounded fraction and melt at final timestep
	hh = ncread(fname, 'h', [1, 1, 1], [Inf,Inf,Inf]); %ice thickness
	hf{i} = squeeze(hh(:,:,end));
	
	ss = ncread(fname, 's', [1, 1, 1], [Inf,Inf,Inf]); %ice surface
	sf{i} = squeeze(ss(:,:,end));

	mm = ncread(fname, 'm', [1, 1, 1], [Inf,Inf,Inf]); %ice melt rate
	mf{i} = squeeze(mm(:,:,end));

	gg = ncread(fname, 'grfrac', [1, 1, 1], [Inf,Inf,Inf]); %ice grounded fraction
	grf{i} = squeeze(gg(:,:,end));

	%comput the bottom
	bb = ss - hh;
	bf{i} = bb;

end 

%% Initialize plots
positions = [0.06, 0.79, 0.64, 0.2; %bathymetry
             0.06, 0.46, 0.64, 0.30; %ice shelf geometry
             0.71, 0.46, 0.1, 0.3;   %forcing temp
             0.82, 0.46, 0.1, 0.3;   %forcing salinity
             0.06, 0.08, 0.35,  0.29;  %zoomed in geometry
             0.415 , 0.08, 0.12, 0.29; %thermal forcing squared
             0.61, 0.08, 0.33, 0.29]; %mean melt rate as a function of M

figure(1); clf;
    
fs = 14; %fontsize
bfs = 16; %big fontsize
icecolor = [194, 227, 236]/255; %ice color
purple = [148,0,211]/255;
blue   = [0, 33, 153]/255;
red    = [153, 0,33]/255;
grey   = [220,220,220]/256;

for i = 1:7
    ax(i) = subplot('Position', positions(i,:));
    hold(ax(i), 'on');
    box(ax(i), 'on')
    ax(i).FontSize = fs;
    ax(i).FontName = 'GillSans';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Make panel a: bed topography %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extend the bed to accomodate ocean region
bedext = zeros(360,50);
bedext(1:300,:) = bed;
bedext(301:end,:) = repmat(bed(end,:), [60,1]);
xext = 500:1e3:(360*1e3 - 500);

% colourmap
cmap = cmocean('ice', 100);
cmap = cmap(1:90,:);

% make plot
imagesc(ax(1), xext/1e3,y/1e3,bedext'); set(gca, 'YDir', 'normal')
plot(ax(1),ax(1).XLim, [0,0], 'w--', 'linewidth', 1.2) % add cross section line
plot(ax(1),[300,300], ax(1).YLim, 'w', 'linewidth', 1.5) %add ice front line

%tidy
%ax(1).XLabel.String = 'x (km)';
ax(1).XLabel.String = '';
ax(1).YLabel.String = 'y (km)';
colormap(ax(1), cmap);
ax(1).YLim = [-24.75,24.75];
ax(1).XLim = [0,359]; %adjust x and y lims so that box is visible
ax(1).XTick = 0:50:350;
ax(1).XTickLabel = {};
ax(1).CLim = [-1200,-500];

%colorbar
c = colorbar(ax(1));
c.Label.String = 'bed depth (m)';
c.Position(1) = c.Position(1)+0.045;
%c.Position(3) = 0.01; %move colorbar out the way
c.FontSize = fs;
c.Label.FontSize = fs;
c.Ticks = -1200:200:-600;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Make panel b: ice cross sections %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add sea level 
plot(ax(2), xext/1e3, zeros(size(xext)), 'linewidth', 1, 'color', 0.5*[1,1,1]);

% add the bed topo 
[nx,ny] = size(h);
fill(ax(2), [xext'/1e3; flip(xext')/1e3], [-1500*ones(size(xext')); flip(bedext(:,floor(ny/2)))], [220,220,220]/256, 'LineStyle','none' ,'FaceAlpha', 0.75);

% add the ice geometry stage one
sp =  smooth(s(:,floor(ny/2))); bp = smooth(b(:,floor(ny/2)));
plot(ax(2),x/1e3,sp, 'k--', 'linewidth', 1.2);
plot(ax(2), x/1e3, bp, 'k--', 'linewidth', 1.2); %pre-melt application
plot(ax(2), x(end)/1e3*[1,1], [bp(end), sp(end)], 'k--', 'linewidth', 1.2)


% add the ice geometry stage two
colmap = parula(sz(2) + 1);
for i = 1:sz(2)
	ss = cell2mat(sf(i));
	bb = cell2mat(bf(i));
	
	%fill in the ice
	if i == sz(2)
		vx = [x/1e3; flip(x)/1e3];
		vy = [bb(:,floor(ny/2)); flip(ss(:,floor(ny/2)))];
		fill(ax(2), vx,vy, icecolor,'linestyle', 'none', 'FaceAlpha', 1)
	end
	
    sp =  smooth(ss(:,floor(ny/2)));
    bp =  smooth(bb(:,floor(ny/2)));
	plot(ax(2), x/1e3,sp, 'color', colmap(i,:), 'linewidth', 1.2);
	plot(ax(2), x/1e3,bp, 'color', colmap(i,:), 'linewidth', 1.2);
    plot(ax(2), x(end)/1e3*[1,1], [bp(end), sp(end)], 'color', colmap(i,:), 'linewidth', 1.2)
end

%add the bed on top
plot(ax(2),xext/1e3,bedext(:,floor(ny/2)), 'k', 'linewidth', 1.2);

% tidy 
ax(2).XLim = [0,359]; %adjust x and y lims so that box is visible
ax(2).YLim = [-1200,500];
ax(2).XTick = ax(1).XTick;
ax(2).XLabel.String = 'x (km)';
ax(2).YLabel.String = 'depth (m)';
ax(2).XLabel.Position(2) =  -1.35e3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Make panel c: temperature and salinity profiles %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zz = -1100:20:0;
ta = nan(size(zz));
ta(zz < -700) = 1.2;
ta(zz > -300) = -1;
idx =( (zz <= -300) & (zz >= -700));
ta(idx) = 1.2 - (2.2)*(zz(idx) + 700)/400;

sa = nan(size(zz));
sa(zz < -700) = 34.6;
sa(zz > -300) = 34.0;
sa(idx) = 34.6 - (0.6)*(zz(idx) + 700)/400;

T0 = 8.32e-2;
lambda = 7.61e-4;
Gamma = 5.73e-2;
Tf = T0 + lambda*zz - Gamma*sa;
dT2 = (ta - Tf).^2;

plot(ax(3), ta, zz, 'color', red, 'linewidth', 2);
plot(ax(3), [-1.2, 1.5], [0,0],'linewidth', 1, 'color', 0.5*[1,1,1]) %add sea level
fill(ax(3), [-1.2, 1.5, 1.5, -1.2], [-1200, -1200, -1100, -1100], grey, 'LineStyle','none' ,'FaceAlpha', 0.75) %add seabed
plot(ax(3),  [-1.2, 1.5], [-1100,-1100], 'k', 'linewidth', 1);

ax(3).YLim = ax(2).YLim;
ax(3).YTick = ax(2).YTick;
ax(3).YTickLabel = {};
ax(3).XLabel.Interpreter = 'latex';
ax(3).XLabel.String = '$T_a$';
ax(3).XLim = [-1.2,1.5];

plot(ax(4), sa, zz, 'color', blue, 'linewidth', 2);
plot(ax(4), [33.8, 34.7], [0,0],'linewidth', 1, 'color', 0.5*[1,1,1]) %add sea level
fill(ax(4), [33.8, 34.7, 34.7, 33.8], [-1200, -1200, -1100, -1100], grey, 'LineStyle','none' ,'FaceAlpha', 0.75) %add seabed
plot(ax(4),  [33.8, 34.7], [-1100,-1100], 'k', 'linewidth', 1);

ax(4).YLim = ax(2).YLim;
ax(4).YTick = ax(2).YTick;
ax(4).YTickLabel = {};
ax(4).XLabel.Interpreter = 'latex';
ax(4).XLabel.String = '$S_a$';
ax(4).XTick = [34,34.6];
ax(4).XLim = [33.8,34.7];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Make panel d: ice cross sections zoomed in %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ns = 8; %smooth nmuber
% add the seabed
fill(ax(5), [x/1e3; flip(x)/1e3], [-1200*ones(size(x)); flip(bed(:,floor(ny/2)))], [220,220,220]/256, 'linestyle', 'none','FaceAlpha', 0.75);

% add the stage 1 calibration ice base
idx = (gr(:, floor(ny/2)) < 1);
sp = x(idx)/1e3; 
bp = smooth(b(idx,floor(ny/2)));
plot(ax(5), sp, bp, 'k--', 'linewidth', 1.2); %pre-melt application
plot(ax(5), x(end)/1e3*[1,1], [bp(end), sp(end)], 'k--', 'linewidth', 1.2)


colmap = parula(sz(2) + 1);
%fill in the ice
i = 1; sz(2);
ss = cell2mat(sf(i));
bb = cell2mat(bf(i));
grfrac = cell2mat(grf(i));

vx = [x/1e3; flip(x)/1e3];
vy = [(bb(:,floor(ny/2))); flip(smooth(ss(:,floor(ny/2))), ns)];
fill(ax(5), vx,vy, icecolor,'linestyle', 'none', 'FaceAlpha', 1)

for i = 1:sz(2)

	ss = cell2mat(sf(i));
	bb = cell2mat(bf(i));
    grfrac = cell2mat(grf(i));

    idx = (grfrac(:, floor(ny/2)) < 1);
    sp =  smooth(ss(idx,floor(ny/2)), ns);
    bp =  smooth(bb(idx,floor(ny/2)), ns);
	plot(ax(5), x(idx)/1e3,sp, 'color', colmap(i,:), 'linewidth', 1.2);
	plot(ax(5), x(idx)/1e3,bp, 'color', colmap(i,:), 'linewidth', 1.2);
    plot(ax(5), x(end)/1e3*[1,1], [bp(end), sp(end)], 'color', colmap(i,:), 'linewidth', 1.2)
	
end

%add the bed
plot(ax(5), x/1e3,bed(:,floor(ny/2)), 'k', 'linewidth', 1);
ax(5).XLabel.String = 'x (km)';
ax(5).YLabel.String = 'depth (m)';
ax(5).XTick = 260:10:300;
ax(5).YTick = -700:200:-300;
ax(5).XLim = [265,305];
ax(5).YLim = [-800,-250];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Make panel e: thermal forcing squared %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(ax(6), dT2, zz, 'color', purple, 'linewidth', 2);
ax(6).YLim = ax(5).YLim;
ax(6).YTick = ax(5).YTick;
ax(6).YTickLabel = {};
ax(6).XLabel.Interpreter = 'latex';
ax(6).XLabel.String = '$(\Delta T)^2$';
ax(6).XLim = [0,15];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Make panel f: mean melt rates %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ms = 0.5:0.25:1.5; %M values
mbar = nan(size(Ms));
for i = 1:sz(2)
	mm = cell2mat(mf(i));
	gg = cell2mat(grf(i));
	mbar(i) = mean(mean(mm(gg ==0))); %mean melt rate in fully floating cells
end

% first plot what melt rate would be if the geometry didn'tchange
yy = mbar(3)*Ms;
plot(ax(7), Ms,yy, 'ko--', 'linewidth',1.5);

plot(ax(7), Ms, mbar, 'k', 'linewidth',1.5);

%add individual points in colours
for i = 1:sz(2)
	plot(ax(7), Ms(i), mbar(i), 'o', 'markersize', 8, 'markeredgecolor', 'k', 'markerfacecolor', colmap(i,:));
end

ax(7).XLabel.Interpreter = 'latex';
ax(7).XLabel.String = '$M$';
ax(7).YLabel.String = 'mean melt rate (m/yr)';
ax(7).XTick = Ms;
ax(7).XTickLabel = {'0.5', '0.75', '1', '1.25', '1.5'};
ax(7).YTick = 10:10:40;
ax(7).YLim  = [10,42];


fig = gcf; 
fig.Position(3:4) = [1080, 520];

%% add an axis for the inset, showing main on different axes
axnew = axes; hold on; box on
axnew.Position = [0.64, 0.26, 0.12, 0.095];

plot(axnew, Ms, mbar, 'k', 'linewidth',1.5);
%add individual points in colours
for i = 1:sz(2)
	plot(axnew, Ms(i), mbar(i), 'o', 'markersize', 6, 'markeredgecolor', 'k', 'markerfacecolor', colmap(i,:));
end


axnew.XLabel.Interpreter = 'latex';
axnew.XLabel.String = '$M$';
axnew.XLabel.String = '';
%axnew.YLabel.String = 'mean melt rate (m/yr)';
axnew.YLabel.String = '';
axnew.XTick = Ms;
axnew.XTickLabel = {'0.5', '0.75', '1', '1.25', '1.5'};
axnew.YTick = [21,24];
axnew.YLim  = [21,24];
axnew.FontName = 'GillSans';
axnew.FontSize = 12;
shg