% Make figure 3 of the attribution manuscript, showing: (a) the bed,  (c) mean
% melt rate as a function of M. (d,e) side view of the steady state
% configuration following calibration and (leave space for a schematic of
% the temperature and salinity profiles)
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

%%
%
% Make the plot
%
figure(1); clf;
positions = [0.1, 0.54, 0.35, 0.35;
             0.585, 0.54, 0.35, 0.35;
             0.1, 0.1, 0.35, 0.35;
              0.585, 0.1, 0.35, 0.35];
%
% (a) Bathymetry
%  

% extend the bed
bedext = zeros(360,50);
bedext(1:300,:) = bed;
bedext(301:end,:) = repmat(bed(end,:), [60,1]);
xext = 500:1e3:(360*1e3 - 500);

ax(1) = subplot('Position', positions(1,:));box on;
cmap = cmocean('ice', 100);
cmap = cmap(1:90,:);
hold on
%contourf(x/1e3,y/1e3,bed', 50, 'linestyle', 'none');
imagesc(xext/1e3,y/1e3,bedext'); set(gca, 'YDir', 'normal')
xl = xlabel('x (km)');
yl = ylabel('y (km)');
colormap(ax(1), cmap);
ax(1).YLim = [-24.75,24.75];
ax(1).XLim = [0,359]; %adjust x and y lims so that box is visible
c = colorbar;
c.Label.String = 'bed depth (m)';

% add cross section line
plot(ax(1).XLim, [0,0], 'w--', 'linewidth', 1.5)

%add ice front lice
plot([300,300], ax(1).YLim, 'w', 'linewidth', 1.5)

% tidy
ax(1).XTick = 0:50:350;
ax(1).CLim = [-1200,-500];
c.Position(1) = 0.455;
c.Position(3) = 0.01;
c.FontSize = 14;
c.Label.FontSize = 16;
c.Ticks = -1200:200:-600;

% 
% (d) Cross sections 
%
[nx,ny] = size(h);
ax(2) = subplot('Position', positions(3,:));box on;
hold on; box on;
%ax(2).YLim = [-24.5,25];
ax(2).XLim = [0,359]; %adjust x and y lims so that box is visible
ax(2).YLim = [-1200,500];
fill([xext'/1e3; flip(xext')/1e3], [-1500*ones(size(xext')); flip(bedext(:,floor(ny/2)))], [220,220,220]/256, 'edgecolor','k' ,'FaceAlpha', 0.75);
ax(2).XTick = 0:50:350;


% add the ice bases
plot(x/1e3, smooth(s(:,floor(ny/2))), 'k--', 'linewidth', 1.2);
plot(x/1e3, smooth(b(:,floor(ny/2))), 'k--', 'linewidth', 1.2); %pre-melt application

colmap = parula(sz(2) + 1);
for i = 1:sz(2)

	ss = cell2mat(sf(i));
	bb = cell2mat(bf(i));
	
	%fill in the ice
	if i == sz(2)
		vx = [x/1e3; flip(x)/1e3];
		vy = [bb(:,floor(ny/2)); flip(ss(:,floor(ny/2)))];
		fill(vx,vy, [0.68,0.85, 0.9],'linestyle', 'none', 'FaceAlpha', 0.75)
	end
	
	plot(x/1e3, smooth(ss(:,floor(ny/2))), 'color', colmap(i,:), 'linewidth', 1.2);
	plot(x/1e3, smooth(bb(:,floor(ny/2))), 'color', colmap(i,:), 'linewidth', 1.2);
	
end

%add the bed
plot(x/1e3,bed(:,floor(ny/2)), 'k');
ax(2).XLabel.String = 'x (km)';
ax(2).YLabel.String = 'depth (m)';




%
% (e) Repeat of (d) with zoomed in on the shelves
% 
ax(3) = subplot('Position', positions(4,:));box on; hold on
fill([x/1e3; flip(x)/1e3], [-1200*ones(size(x)); flip(bed(:,floor(ny/2)))], [220,220,220]/256, 'linestyle', 'none','FaceAlpha', 0.75);


% add the ice bases
idx = (gr(:, floor(ny/2)) < 1);
plot(x(idx)/1e3, smooth(b(idx,floor(ny/2))), 'k--', 'linewidth', 1.2); %pre-melt application

colmap = parula(sz(2) + 1);
%fill in the ice
i = sz(2);
ss = cell2mat(sf(i));
bb = cell2mat(bf(i));
grfrac = cell2mat(grf(i));

vx = [x/1e3; flip(x)/1e3];
vy = [(bb(:,floor(ny/2))); flip(smooth(ss(:,floor(ny/2))))];
fill(vx,vy, [0.68,0.85, 0.9],'linestyle', 'none', 'FaceAlpha', 0.75)

for i = 1:sz(2)

	ss = cell2mat(sf(i));
	bb = cell2mat(bf(i));
    grfrac = cell2mat(grf(i));
	
	
	idx = grfrac(:,floor(ny/2)) < 1;
	plot(x(idx)/1e3, smooth(ss(idx,floor(ny/2))), 'color', colmap(i,:), 'linewidth', 1.4);
	plot(x(idx)/1e3, smooth(bb(idx,floor(ny/2)),8), 'color', colmap(i,:), 'linewidth', 1.4);
	
end

%add the bed
plot(x/1e3,bed(:,floor(ny/2)), 'k');
ax(3).XLabel.String = 'x (km)';
ax(3).YLabel.String = 'depth (m)';
ax(3).XTick = 260:10:300;
ax(3).XLim = [265,305];
ax(3).YLim = [-750,-350];


%
% (b) Mean shelf melt rate
%
ax(4) = subplot('Position', positions(2,:));box on;
hold on; box on;
gammas = [8,9,10,11,12]*1e-3;
mbar = nan(size(gammas));
for i = 1:sz(2)
	mm = cell2mat(mf(i));
	gg = cell2mat(grf(i));
	mbar(i) = mean(mean(mm(gg ==0))); %mean melt rate in fully floating cells
end
plot(gammas*1e3, mbar, 'k--', 'linewidth',1.5);

%add individual points in colours
for i = 1:sz(2)
	plot(gammas(i)*1e3, mbar(i), 'o', 'markersize', 8, 'markeredgecolor', 'k', 'markerfacecolor', colmap(i,:));
end
ax(4).XLabel.Interpreter = 'latex';
ax(4).XLabel.String = '$M$';
ax(4).YLabel.String = 'mean melt rate (m/yr)';
ax(4).XTick = 8:12;
ax(4).XTickLabel = {'0.5', '0.75', '1', '1.25', '1.5'};
ax(4).YTick = 20:25;
ax(4).YLim  = [21,24];
ax(4).Position(1) = 0.75; %shift over to make room for schematics of forcing
ax(4).Position(3) = ax(3).Position(1) + ax(3).Position(3) - ax(4).Position(1); %make right hand sides align


for i = 1:4
    ax(i).FontSize = 14;
    ax(i).FontName = 'GillSans';
end


fig = gcf; 
fig.Position(3:4) = [1080, 520];
