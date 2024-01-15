% Make figure 2 of the short attribution manuscript, showing: 
% (a) the bed elevation
% (b) the ice geometry in initial condition
% (c) forcing profiles
% (d) forcing of a single member
% (e) sea level rise as a function of time for this realization and different values of M
% (f) forcing profiles of all members
% (g) sea level rise at t = 100 as a function of M and different ensembles

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

end 

%% Initialize plots
positions = [0.06, 0.82, 0.64, 0.17; %bathymetry
             0.06, 0.61, 0.64, 0.19; %ice shelf geometry
             0.71, 0.61, 0.1,  0.19; %forcing temp
             0.82, 0.61, 0.1,  0.19; %forcing salinity
             0.06, 0.36, 0.39, 0.18; %forcing relization
             0.06, 0.07, 0.39, 0.23; %slr for this different forcing
             0.53, 0.36, 0.39, 0.18; %all realizations of forcing
             0.53, 0.07, 0.39, 0.23]; % final slr at different times and M and ensembles

figure(1); clf;
    
fs = 14; %fontsize
bfs = 16; %big fontsize
icecolor = [194, 227, 236]/255; %ice color
purple = [148,0,211]/255;
blue   = [0, 33, 153]/255;
red    = [153, 0,33]/255;
grey   = [220,220,220]/256;

%colourmap for anthro vs counterfactual
colmapg = nan(2,3); 
colmapg(1,:) = [255,152,51]/255;  %anthro
colmapg(2,:) = [0,153, 153]/255; %counter
% for poster
% colmapg(2,:) = [7,54,125]/255; %dark blue
% colmapg(1,:) = [248,200,44]/255; %yellow


for i = 1:8
    ax(i) = subplot('Position', positions(i,:));
    hold(ax(i), 'on');
    box(ax(i), 'on')
    ax(i).FontSize = fs;
    ax(i).FontName = 'Arial';
end
%%
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

% % add the ice geometry stage one
% sp =  smooth(s(:,floor(ny/2))); bp = smooth(b(:,floor(ny/2)));
% plot(ax(2),x/1e3,sp, 'k--', 'linewidth', 1.2);
% plot(ax(2), x/1e3, bp, 'k--', 'linewidth', 1.2); %pre-melt application
% plot(ax(2), x(end)/1e3*[1,1], [bp(end), sp(end)], 'k--', 'linewidth', 1.2)


% add the ice geometry stage two
colmap = parula(sz(2) + 1);
colmap = zeros(sz(2)+1,3); %make all black
for i = sz(2)
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
%%%%%%%%%%%% Make panel c,d: temperature and salinity profiles %%%%%%%%%%%%
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
%%%%%%%%%%%% Make panel e: single value of forcing %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ensemble      = 1; %1: anthro trend, 2: no trend
member        = 5; %.e.g.

% load in the data
fa = load('data/forcing_anomalies.mat');
pc = fa.ss(ensemble,member).pc;
tt  = fa.ss(ensemble,member).t;

% sub-sample and make positive
pcs = pc(1:1:end);
ts  = tt(1:1:end);
pcspos = pcs;
pscpos(pcs < -500) = -500;
pcsneg = pcs;
pcsneg(pcs > -500) = -500;


%fill positive
xf =  [ts,flip(ts)];
yf =  [-500*ones(size(pcspos)); flip(pcspos)];
fill(ax(5), xf, yf, [153, 33, 33]/255, 'LineStyle', 'none', 'FaceAlpha',0.7);

%fill negative
yf =  [pcsneg;  -500*ones(size(pcspos)), ];
fill(ax(5), xf, yf, [0, 33, 153]/255, 'LineStyle', 'none', 'FaceAlpha',0.7);

%plot(ax(5), ts,  pcs, 'linewidth', 1.4, 'color', 'k');

%plot counterfactual trend
plot(ax(5), [0,100] ,-500*[1,1],'k', 'linewidth', 1, 'Color',  'k')

%plot anthro trend
%plot(ax(5), [0,100] ,[-500, -400],'k', 'linewidth', 2)


%tidy 
ax(5).XLim = [0,100];
ax(5).XTick = 0:20:100;
%ax(5).XTickLabel = {};
ax(5).YLim = [-800, -200];
ax(5).YTick = [-700,-500,-300];
ax(5).YLabel.String = '$P_c(t)$';
ax(5).XLabel.String = 'time (years)'; 
ax(5).YLabel.Interpreter = 'latex';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Make panel f: SLR for the forcing in e %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% preliminary stuff
ice_data = load('data/WAVI-ensemble-data.mat'); ice_data = ice_data.ss; %load the ice sheet run data in
dx = 1e3;
dy = 1e3;
rhoi = 918;  %ice density
rhow = 1028; %water density
fpath = strcat('data/ATTR_00000_outfile.nc');
bed   = ncread(fpath, 'b', [1, 1, 1], [Inf, Inf, 1]);
float_thick = abs(rhow/rhoi *bed); %thickness at which floatation occurs
%get slr data

%
% get slr curve for each M
%
Ms        = 1:5; %correspond to M = 0.5:0.25:1.5
Ms_act    = 0.5:0.25:1.5;
colmap = (cmocean('ice',2*length(Ms)+10)); %coloumap
colmap = colmap(4:2:end-6,:);
colmap = flipud(colmap);

for iM = 1:length(Ms)
    %get slr at associataed with these
    hh = ice_data(iM,ensemble, member).h; %ice thickness
    idx = hh > float_thick;
    dh = hh - float_thick;
    dh(~idx) = 0;               %mask anything below float thick to 0
    vaf = squeeze(sum(sum(dh,2),1)*dx*dy); %vaf as a function of time
    slr = (vaf(1) - vaf)/ 395 / 1e9; %SLR in mm at all times 
    tt = ice_data(iM,ensemble,member).t;
    ss = smooth(slr,5 );
    p = plot(ax(6), tt, ss, 'linewidth', 2, 'Color', colmap(iM,:));
    
    %add a final point
    scatter(ax(6), tt(end), ss(end),60,colmap(iM,:), 'Filled', 'MarkerEdgeColor', 'none' )
    
    finval(iM) = ss(end);
end



ax(6).XLim = [0,100];
ax(6).XLabel.String = 'time (years)';
ax(6).YLabel.String = '\Delta SLR (mm)';
ax(6).YLim = [-1,4];

%make agree with e
ax(6).XTick = ax(5).XTick;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Make panel g: all realizations of forcing %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('data/forcing_anomalies.mat');

ax(1) = subplot('Position', positions(1,:)); hold on; box on
pc_mean = zeros(2001,2);
sz = [2,30]; %size of the ensemble
for ie = 1:2
    for im = 1:sz(2)
        tt = ss(ie,im).t;
        pp = ss(ie, im).pc;
        p = plot(ax(7), tt(1:20:end), pp(1:20:end), 'linewidth', 1, 'HandleVisibility','off');
        p.Color = [colmapg(ie,:),0.15];
        pc_mean(:,ie) = pc_mean(:,ie) + pp;
        
    end

end

pc_mean = pc_mean/sz(2);
for ie = 1:2   
    plot(ax(7), tt(1:20:end),pc_mean(1:20:end,ie), '--',  'linewidth', 1.75, 'color', colmapg(ie,:));
    % add the running trend 
    plot(ax(7), tt(1:20:end),smooth(pc_mean(1:20:end,ie),10), '--', 'linewidth', 1.75, 'color', colmapg(ie,:));
end
%plot(ax(7), tt(1:20:end),-500 + tt(1:20:end), '--','linewidth', 1.25, 'color', colmapg(1,:), 'HandleVisibility','off'); %add the anthro trend
%plot(ax(7), tt(1:20:end),-500*ones(size(tt(1:20:end))), '--','linewidth', 1.25, 'color', colmapg(2,:), 'HandleVisibility','off'); %add the natural trend

ax(7).XLim = [0, 100];
ax(7).YLim = [-700, -300];
ax(7).XLabel.String = 'time (yrs)';
ax(7).YLabel.String = 'pycnocline depth (m)';
%ylim([-750, -250]);
ax(7).YTick = -900:200:-100;
ax(7).XTick = 0:20:100;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Make panel h: SLR at t = 100 for different M and realization %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ms        = 1:5; %correspond to M = 0.5:0.25:1.5
Ms_act    = 0.5:0.25:1.5;
ensembles = 1:2; %1: anthro ensemble, 2: counterfactual ensemble
members   = 1:6; 
t_out     = 100;  %output time

slrs      = nan(length(ensembles), length(members), length(Ms));
for ie = 1:length(ensembles)
for im = 1:length(members)
    for iM = 1:length(Ms)
        %get slr at associataed with these
        hh = ice_data(iM,ensembles(ie),im).h; %ice thickness
        idx = hh > float_thick;
        dh = hh - float_thick;
        dh(~idx) = 0;               %mask anything below float thick to 0
        vaf = squeeze(sum(sum(dh,2),1)*dx*dy); %vaf as a function of time
        slr = (vaf(1) - vaf)/ 395 / 1e9; %SLR in mm at all times 

        %reduce to single time output
        tt = ice_data(iM,ensemble,im).t;
        [~,idx] = min(abs(tt - t_out));
        slrs(ie, im,iM) = slr(idx);
        

    end
end
end

lw = 1.5;

%add cointerfactual results
for i = 1:length(members)
     
    plot(ax(8), Ms_act, squeeze(slrs(2,i,:)), '-','color',  colmapg(2,:), 'markerfacecolor',  colmap(2,:), 'linewidth', lw)
end

%add anthro results
for i = 1:length(members)
    plot(ax(8), Ms_act, squeeze(slrs(1,i,:)), '-','color', colmapg(1,:), 'markerfacecolor', colmap(1,:), 'linewidth', lw)
end



ax(8).XTick = Ms_act;
ax(8).XLabel.Interpreter = 'latex';
ax(8).XLabel.String = '$M$';
ax(8).YTick = 0:4;
ax(8).YLabel.String = 'sea level rise (mm)';
ax(8).YLim = [-0.1, 4];

% add the points in (g)
for iM = 1:length(Ms)
    scatter(ax(8),Ms_act(iM), finval(iM),80,colmap(iM,:), 'Filled' , 'MarkerEdgeColor', 'none')
end
%legend(ax(8), {'anthropogenic trend', 'no trend'}, 'FontSize', 16, 'Location', 'SouthEast')



%% final tidying

fig = gcf; 
fig.Position(3:4) = [1080, 710];
