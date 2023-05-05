% Make supplementary figure "C" of the attribution manuscript, showing: 
% (a) schematic evolution of grounded volume during the calibration
% procedure [NB: this is produced externally]
% (b) mean melt rate after second calibration phase as a function of M
% (c) ice shelf geometry after calibration procedure
% (d) squared ice shelf forcing.

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
dx = 1e3;
dy = 1e3; 

%% Plot setup
positions = [0.10, 0.60, 0.40, 0.39; %VAF evolution 
             0.55, 0.60, 0.40, 0.39;
             0.10, 0.10, 0.64, 0.39; %ice shelf geometric
             0.75, 0.10, 0.20, 0.39]; %forcing profiles

figure(1); clf;
    
fs = 14; %fontsize
icecolor = [194, 227, 236]/255; %ice color
purple = [148,0,211]/255;
blue   = [0, 33, 153]/255;
red    = [153, 0,33]/255;
grey   = [220,220,220]/256;
Mcolmap = parula(6); %different shelf colours

for i = 1:4
    ax(i) = subplot('Position', positions(i,:));
    hold(ax(i), 'on');
    box(ax(i), 'on')
    ax(i).FontSize = fs;
    ax(i).FontName = 'Arial';
end
%% Fetch the data
%
% Get the results for the no-melt stage of the calibration
%
fpath = strcat('data/ATTR_00000_outfile.nc');
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
tf       = cell(sz);  %times
grv      = cell(sz);  %grounded volume values
for i = 1:sz(2)
	%filename
	fname =  strcat('data/ATTR_', run_nums(i), '_outfile.nc');

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
    
    tf{i} = ncread(fname, 'TIME');

    grv{i} =  squeeze(sum(sum(hh.*gg,1),2)*dx*dy);

end 

%% Make panel c: ice cross sections zoomed in on shelf

ns = 8; %smooth nmuber
ny = 50; %number of grid cells in y
% add the seabed
fill(ax(3), [x/1e3; flip(x)/1e3], [-1200*ones(size(x)); flip(bed(:,floor(ny/2)))], [220,220,220]/256, 'linestyle', 'none','FaceAlpha', 0.75);

% add the stage 1 calibration ice base
idx = (gr(:, floor(ny/2)) < 1);
sp = x(idx)/1e3; 
bp = smooth(b(idx,floor(ny/2)));
plot(ax(3), sp, bp, 'k--', 'linewidth', 1.2); %pre-melt application
plot(ax(3), x(end)/1e3*[1,1], [bp(end), sp(end)], 'k--', 'linewidth', 1.2)


colmap = parula(sz(2) + 1);
%fill in the ice
i = 1; sz(2);
ss = cell2mat(sf(i));
bb = cell2mat(bf(i));
grfrac = cell2mat(grf(i));

vx = [x/1e3; flip(x)/1e3];
vy = [(bb(:,floor(ny/2))); flip(smooth(ss(:,floor(ny/2))), ns)];
fill(ax(3), vx,vy, icecolor,'linestyle', 'none', 'FaceAlpha', 1)

for i = 1:sz(2)

	ss = cell2mat(sf(i));
	bb = cell2mat(bf(i));
    grfrac = cell2mat(grf(i));

    idx = (grfrac(:, floor(ny/2)) < 1);
    sp =  smooth(ss(idx,floor(ny/2)), ns);
    bp =  smooth(bb(idx,floor(ny/2)), ns);
	plot(ax(3), x(idx)/1e3,sp, 'color', colmap(i,:), 'linewidth', 1.5);
	plot(ax(3), x(idx)/1e3,bp, 'color', colmap(i,:), 'linewidth', 1.5);
    plot(ax(3), x(end)/1e3*[1,1], [bp(end), sp(end)], 'color', colmap(i,:), 'linewidth', 1.5)
	
end

%add the bed
plot(ax(3), x/1e3,bed(:,floor(ny/2)), 'k', 'linewidth', 1);
ax(3).XLabel.String = 'x (km)';
ax(3).YLabel.String = 'depth (m)';
ax(3).XTick = 260:10:300;
ax(3).YTick = -700:200:-300;
ax(3).XLim = [265,305];
ax(3).YLim = [-800,-250];

%% Make panel d: thermal forcing squared
%
% generate profiles
%
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

%
% make the plot
%
plot(ax(4), dT2, zz, 'color', purple, 'linewidth', 2);
ax(4).YLim = ax(3).YLim;
ax(4).YTick = ax(3).YTick;
ax(4).YTickLabel = {};
ax(4).XLabel.Interpreter = 'latex';
ax(4).XLabel.String = '$(\Delta T)^2$';
ax(4).XLim = [0,15];

%% Make panel b: mean melt rates

Ms = 0.5:0.25:1.5; %M values
mbar = nan(size(Ms));
for i = 1:sz(2)
	mm = cell2mat(mf(i));
	gg = cell2mat(grf(i));
	mbar(i) = mean(mean(mm(gg ==0))); %mean melt rate in fully floating cells
end

% first plot what melt rate would be if the geometry didn'tchange
yy = mbar(3)*Ms;
plot(ax(2), Ms,yy, 'ko--', 'linewidth',1.5);

plot(ax(2), Ms, mbar, 'k', 'linewidth',1.5);

%add individual points in colours
for i = 1:sz(2)
	plot(ax(2), Ms(i), mbar(i), 'o', 'markersize', 8, 'markeredgecolor', 'k', 'markerfacecolor', colmap(i,:));
end

ax(2).XLabel.Interpreter = 'latex';
ax(2).XLabel.String = '$M$';
ax(2).YLabel.String = 'mean melt rate (m/yr)';
ax(2).XTick = Ms;
ax(2).XTickLabel = {'0.5', '0.75', '1', '1.25', '1.5'};
ax(2).YTick = 10:10:40;
ax(2).YLim  = [5,35];


fig = gcf; 
fig.Position(3:4) = [1080, 520];

%% add an axis for the inset, showing main on different axes
axnew = axes; hold on; box on
axnew.Position = [0.76, 0.66, 0.18, 0.13];

%add mitgcm mean melt rate
plot([0.5, 1.5], [23,23], 'r', 'linewidth', 2);

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
axnew.FontName = 'Arial';
axnew.FontSize = 12;
shg