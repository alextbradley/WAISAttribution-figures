% Make figure 3 of the attribution manuscript
%
% Stochastic retreat: grounded volume as a function of time for gamma_T = 10e-3 with natural and anthro and corresponding forcing profiles
%
%
% Preliminaries
%
addpath('plottools')
gendata = 1; %gendata loop flag

%
% Constants
%
rhoi = 918;  %ice density
rhow = 1028; %water density
dx   = 1000;
dy   = 1000; %grid resolution

%
% setup plot
%
positions = [0.08, 0.15, 0.4, 0.8;
    0.55, 0.15, 0.4, 0.8];
%
% Part 1: plot of stochastic forcing
%
load('data/forcing_anomalies.mat');
figure(1); clf;
colmap = [0, 33, 153; 153,0,63]/255;
colmap = lines(2); %colormap to differentiate between anthro and non
colmap(1,:) = [ 0,    0.45    0.84];
colmap(2,:) = [ 0.80    0.24    0.1];

ax(1) = subplot('Position', positions(1,:)); hold on; box on
pc_mean = zeros(2001,2);
sz = [2,20]; %size of the ensemble
for ie = 1:sz(1)
    for im = 1:sz(2)
        tt = ss(ie,im).t;
        pp = ss(ie, im).pc;
        p = plot(tt, pp, 'linewidth', 1, 'HandleVisibility','off');
        p.Color = [colmap(ie,:),0.15];
        pc_mean(:,ie) = pc_mean(:,ie) + pp;

    end
    pc_mean = pc_mean/sz(2);
    plot(tt,pc_mean(:,ie), 'linewidth', 1.75, 'color', colmap(ie,:));
end
plot(tt,-500 + tt, '--','linewidth', 1.25, 'color', colmap(1,:), 'HandleVisibility','off'); %add the anthro trend
plot(tt,-500*ones(size(tt)), '--','linewidth', 1.25, 'color', colmap(2,:), 'HandleVisibility','off'); %add the natural trend

xlim([0, 100]);
xlabel('time (yrs)');
ylabel('pycnocline depth (m)');
%ylim([-750, -250]);
ax(1) = gca; ax(1).FontSize = 11;
ax(1).YTick = -900:200:-100;
legend({'anthropogenic trend', 'no trend'}, 'FontSize', 16, 'Location', 'SouthEast')

%%
% Part 2: sea level rise

%
% Run info
%
gammas    = [11; 11];
gidx      = 4;          %index in the values: (4 corresponds to 11e-4)
ensembles = [1; 2]; %anthro (with trend) first, then natural (no trend)
members   = 1:20;
members   = repmat(members, [2,1]);
sz        = size(members);

%
% Initialize storage
%
if gendata
    ss = load('data/WAVI-ensemble-data.mat');
    s = ss.ss;
    vaf  = cell(sz);
    t    = cell(sz); %times corresponding to vaf (different to those for forcing)
    dx   = 1000;
    dy   = 1000;

    %
    % compute the ice bed
    %
    fpath = strcat('data/ATTR_00000/outfile.nc');
    bed   = ncread(fpath, 'b', [1, 1, 1], [Inf, Inf, 1]); %bed topo
    float_thick = abs(rhow/rhoi *bed);

    %
    % Get the results for VAF and associate time and forcing and associated time
    %
    for i = 1:sz(1) %forcing type
        for j = 1:sz(2) %member number

        	hh = s(gidx,i,j).h; %ice thickness
            idx = hh > float_thick;
            dh = hh - float_thick;
            dh(~idx) = 0;
            vv = sum(sum(dh,2),1)*dx*dy;
            
%           gg = s(i,j).grfrac; %ice grounded fraction
%         
%          	vv = hh.*gg;
%         	vv = sum(sum(vv,2),1)*dx*dy;
%         	vv = squeeze(vv);
        	vaf{i,j} = vv;
            tt = s(i,j).t;
            t{i,j} = tt;

        end %end loop over member number
    end %end loop over ensembles
end %end gendata flag

%
% Make plot
%
ax(2) = subplot('Position', positions(2,:)); hold on; box on
for i = 1:sz(1)
    for j = 1:sz(2)
    	tt = cell2mat(t(i,j));
    	vv = cell2mat(vaf(i,j));
    	volume_change = vv(1) - vv;
    	SLR = volume_change / 395 / 1e9; %SLR in mm
    	p = plot(tt,smooth(SLR, 4),'linewidth', 1.25);


    	p.Color = [colmap(i,:),0.6];
    end
end


xlim([0, 100]);
xlabel('time (yrs)');
ylabel('sea level rise (mm)');
ax(2) = gca;
ax(2).FontSize = ax(1).FontSize;
ax(2).YLim = [-1,4];
grid on

fig = gcf;
fig.Position(3:4) = [950, 360];

%%
for i = 1:2
    ax(i).FontSize = 16;
    ax(i).FontName = 'GillSans';
end
fig = gcf;
fig.Position(3:4) = [1080, 380];
