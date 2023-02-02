% Make figure 2 of the attribution manuscript, showing:
% (a) the bed topography
% (b) side view of the steady state configuration
% (c) mean melt rate as a function of \gamma.
%
% Alex Bradley, 10/11/22. MIT License.

%
% Preliminaries
%
addpath('plottools')
outfile_path = '/data/icesheet_output/aleey/wavi/'; %change to full path of result location

%
% Constants
%
rhoi = 918;  %ice density
rhow = 1028; %water density

%
% Get the results for the no-melt stage of the calibration
%
fpath = strcat(outfile_path, 'ATTR_00000/run/outfile.nc');
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
bf 	 = cell(sz);  %ice bottom
for i = 1:sz(2)
	%filename
	fname =  strcat(outfile_path, 'ATTR_', run_nums(i), '/run/outfile.nc');

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


%
% Make the plot
%
figure(1); clf;

%
% (a) Bathymetry
%  
ax(1) = subplot(2,2,1);box on;
cmap = cmocean('haline');
hold on
contourf(x/1e3,y/1e3,bed', 50, 'linestyle', 'none');
%xl = xlabel('x (km)');
yl = ylabel('y (km)');
colormap(ax(1), cmap);
ax(1).YLim = [-25,25];
ax(1).XLim = [-1,300]; %adjust x and y lims so that box is visible
c = colorbar;
c.Label.String = 'bed depth (m)';

% add cross section line
plot(ax(1).XLim, [0,0], 'w--', 'linewidth', 1.5)

% 
% (b) Cross sections 
%
[nx,ny] = size(h);
ax(2) = subplot(2,2,3); hold on; box on;
%ax(2).YLim = [-24.5,25];
ax(2).XLim = [-1,300]; %adjust x and y lims so that box is visible
ax(2).YLim = [-1200,500];
ax(2).Position(3:4) = ax(1).Position(3:4);
fill([x/1e3; flip(x)/1e3], [-1200*ones(size(x)); flip(bed(:,floor(ny/2)))], [220,220,220]/256, 'linestyle', 'none','FaceAlpha', 0.75);


% add the ice bases
colmap = parula(sz(2) + 1);
for i = 1:sz(2)

	ss = cell2mat(sf(i));
	bb = cell2mat(bf(i));
	
	%fill in the ice
	if i == sz(2);
		vx = [x/1e3; flip(x)/1e3];
		vy = [bb(:,floor(ny/2)); flip(ss(:,floor(ny/2)))];
		fill(vx,vy, [0.68,0.85, 0.9],'linestyle', 'none', 'FaceAlpha', 0.75)
	end
	
	plot(x/1e3, ss(:,floor(ny/2)), 'color', colmap(i,:));
	plot(x/1e3, bb(:,floor(ny/2)), 'color', colmap(i,:));
	
end

%add the bed
plot(x/1e3,bed(:,floor(ny/2)), 'k');
ax(2).XLabel.String = 'x (km)';
ax(2).YLabel.String = 'depth (m)';

% fix (a) to be the same
ax(1).XTick = ax(2).XTick;
ax(1).XTickLabel = {};
c.Position(1) = 0.47; %just put the cbar out the way for now

%
% (c) Repeat of (b) with zoomed in on the shelves
% 
ax(3) = subplot(2,2,4); hold on; box on;
ax(3).XLim = [260,299]; %adjust x and y lims so that box is visible
ax(3).YLim = [-1200,500];
ax(3).Position(3:4) = ax(1).Position(3:4);
fill([x/1e3; flip(x)/1e3], [-1200*ones(size(x)); flip(bed(:,floor(ny/2)))], [220,220,220]/256, 'linestyle', 'none','FaceAlpha', 0.75);


% add the ice bases
colmap = parula(sz(2) + 1);
for i = 1:sz(2)

	ss = cell2mat(sf(i));
	bb = cell2mat(bf(i));
	
	%fill in the ice
	if i == sz(2);
		vx = [x/1e3; flip(x)/1e3];
		vy = [bb(:,floor(ny/2)); flip(ss(:,floor(ny/2)))];
		fill(vx,vy, [0.68,0.85, 0.9],'linestyle', 'none', 'FaceAlpha', 0.75)
	end
	
	plot(x/1e3, ss(:,floor(ny/2)), 'color', colmap(i,:));
	plot(x/1e3, bb(:,floor(ny/2)), 'color', colmap(i,:));
	
end

%add the bed
plot(x/1e3,bed(:,floor(ny/2)), 'k');
ax(3).XLabel.String = 'x (km)';
ax(3).YLabel.String = 'depth (m)';

%
% (d) Mean shelf melt rate
%
ax(4) = subplot(2,2,2); hold on; box on;
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
ax(4).XLabel.String = '$\gamma_T \times 10^3$';
ax(4).YLabel.String = 'mean melt rate (m/yr)';
ax(4).Position(1) = 0.605;
ax(4).Position(3) = (ax(3).Position(1) + ax(3).Position(3)) - ax(4).Position(1);
ax(4).XTick = 8:12;
ax(4).YTick = 20:25;
ax(4).YLim  = [20,25];


