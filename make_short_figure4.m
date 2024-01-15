% Make figure 4 of the manuscript, showing (a) the anthropogenic
% enhancement ratio and (b)--(d) example curves along this trajectory
%
% 27/02/23, ATB (aleey@bas.ac.uk), MIT licence
%

%% Preliminaries
% load in the data
pdata = load('data/shortfigure-3data.mat'); %distribution data
wavdat = load('data/WAVI-ensemble-data.mat'); wavdat = wavdat.ss; %for use in the SLR curves
load('data/shortfigure4-bootstrapdata.mat'); % load in bootstrapping data

addpath('plottools');

%% Plot setup
fig = figure(1); clf;
fig.Color = 'w';
fig.Position(3:4) = [1200, 420];

ht = 0.85; %total height
gapy = 0.05;
hts = (ht - 2*gapy)/3;
starty = 0.11;

positions = [0.11, starty, 0.35, 0.85;
    0.57, (starty + 2*hts +2*gapy),  0.35,  hts;
    0.57, (starty + hts +gapy), 0.35, hts;
    0.57, starty, 0.35, hts];
for i = 1:4
    ax(i) = subplot('Position', positions(i,:));
    hold(ax(i), 'on');
    box(ax(i), 'on');
    ax(i).FontName = 'Arial';
    ax(i).FontSize = 14;
end

%
% setup colors
%
colmap = nan(2,3);
colmap(1,:) = [255,152,51]/255;  %anthro
colmap(2,:) = [0,153, 153]/255; %counter

% % for poster
% colmap(2,:) = [7,54,125]/255; %dark blue
% colmap(1,:) = [248,200,44]/255; %yellow

fillcol = [33,0,152]/255; %for scenario examples
tailcol = [255,133,133]/255;

%create colormap between these
T = [colmap(2,:);
    1,1,1;
    colmap(1,:)];

x = [0
    50
    100]; %intervals of colormap (choose middle number to match clims in a)

cmap = interp1(x/100,T,linspace(0,1,255));

% colourmap for significance levels
sigcmap = [ 0.9290    0.6940    0.1250;
    0, 200/255, 197/255];
sigcmap = lines(2);
sigcmap = [231, 76, 60;
    131, 52, 131 ]/255;
%% (a) Anthropogenic enhancement and significance contours
% subsample and smoothing parameters
nxs = 10; %how finely to subsample x
nsmooth = nxs*3;
t            = pdata.t;

%extract the distributions
anth_all           = pdata.vals(:,1,:,:);
nat_all            = pdata.vals(:,2,:,:);
anth_mean          = squeeze(mean(anth_all,3)); %take mean over ensemble members
nat_mean           = squeeze(mean(nat_all,3)); %take mean over ensemble members

%subsample stuff
anth_all_subsamp  = squeeze(anth_all(:,:,:,1:nxs:end));
nat_all_subsamp   = squeeze(nat_all(:,:,:,1:nxs:end));
anth_mean_subsamp = anth_mean(:,1:nxs:end);
nat_mean_subsamp  = nat_mean(:,1:nxs:end);

%smooth stuff
anth_all_subsamp_smooth  = nan(size(anth_all_subsamp));
nat_all_subsamp_smooth   = nan(size(nat_all_subsamp));
anth_mean_subsamp_smooth = nan(size(anth_mean_subsamp));
nat_mean_subsamp_smooth = nan(size(nat_mean_subsamp));


for it = 1:length(t)
    for im = 1:40
        anth_all_subsamp_smooth(it,im,:) = smooth(squeeze(anth_all_subsamp(it,im,:)), nsmooth)';
        nat_all_subsamp_smooth(it,im,:) = smooth(squeeze(nat_all_subsamp(it,im,:)), nsmooth)';
    end

    anth_mean_subsamp_smooth(it,:) = smooth(squeeze(anth_mean_subsamp(it,:)), nsmooth)';
    nat_mean_subsamp_smooth(it,:) = smooth(squeeze(nat_mean_subsamp(it,:)), nsmooth)';

end

% issig_pt01         = pdata.is_significant_pt01(:,1:nxs:end);
% issig_pt1          = pdata.is_significant_pt1(:,1:nxs:end);
% issig_pt05         = pdata.is_significant_pt05(:,1:nxs:end);
xx                 = pdata.x(1:nxs:end);

AER = nan(length(t), length(xx));
AER_smooth = nan(length(t), length(xx));
AER_upper = nan(length(t), length(xx));
AER_lower = nan(length(t), length(xx));
AER_lower_20 = nan(length(t), length(xx));
for i = 1:length(t)

    % smooth distributions
    %store
    AER_smooth(i,:) = anth_mean_subsamp_smooth(i,:)./nat_mean_subsamp_smooth(i,:);
    AER(i,:)        =  anth_mean_subsamp(i,:)./nat_mean_subsamp(i,:);

    AER_upper(i,:) = anth_ci_upper(i,:)./nat_ci_lower(i,:);
    AER_lower(i,:) = anth_ci_lower(i,:)./nat_ci_upper(i,:);

end

%AER(isinf(AER)) = 1e4;
AER_smooth = smooth2a(AER_smooth, 2,10);
p = imagesc(ax(1), xx, t, log10(AER_smooth));

set(p, 'AlphaData', ~isnan(AER));
set(gca, 'YDir', 'normal');
c = colorbar(ax(1));
c.Ticks = -1:0.5:1.5;
c.TickLabels = {'10^{-1}','10^{-0.5}' '10^{0}', '10^{0.5}','10^{1}', '>10^{1.5}'};
colormap(ax(1), cmap)
ax(1).CLim = [-1,1];
ax(1).XLim = [-0.3, 4];
ax(1).YLim = [10,100];
ax(1).YLabel.String = 'time (years)';
ax(1).XLabel.String = 'SLR (mm)';
ax(1).YTick = 20:20:100;

% do a different colour for the tail?
% axnew = axes;
% AERinf = AER;
% AERinf(~isinf(AER)) = nan;
% AERinf(isinf(AER)) = 1;
% AERinf(:,1:230) = nan;  %remove anything at negative slr
%
% p = imagesc(axnew, xx, t, AERinf);
% set(p, 'AlphaData', ~isnan(AERinf))
% set(axnew, 'YDir', 'normal');
% axnew.Visible = 'off';
% colormap(axnew, tailcol)
%
% axnew.XLim = ax(1).XLim;
% axnew.YLim = ax(1).YLim;
% axnew.Position = ax(1).Position;
% linkaxes([ax(1), axnew])

c.Position(1) = 0.47; %move after linking

%
% add significance contorurs
%
AERrad = (AER_lower > 1);
AERrad_20 = (AER_lower_20 > 1);

%contour(ax(1), xx,t,smooth2a(AERrad, 2, nsmooth), [0.5,0.5],'color','k', 'linewidth', 1.5)
%% Make (b) showing different trajectories
%
% trajectory info
%
ims = [9,11,5]; %5as high, 9 or 16 as low, 10 as mid
ims =[19, 23,36];
iM = 5;
ie = 1;

%
% get stuff needed for the calculation
%
fpath = strcat('data/ATTR_00000_outfile.nc');
bed   = ncread(fpath, 'b', [1, 1, 1], [Inf, Inf, 1]); %bed topo
float_thick = abs(1028/918 *bed);
dx = 1e3; dy = 1e3;

%
% ice sheet data
%

for ii = 1:3
    im = ims(ii);
    hh = wavdat(iM,ie,im).h; %ice thickness
    idx = hh > float_thick;
    dh = hh - float_thick;
    dh(~idx) = 0;
    vv = sum(sum(dh,2),1)*dx*dy;
    vv = squeeze(vv);
    volume_change = vv(1) - vv;
    SLR = volume_change / 395 / 1e9; %SLR in mm
    SLR = smooth(SLR);
    tt = wavdat(iM,ie,im).t;
    plot(ax(1), SLR, tt, 'k--', 'linewidth', 1.75, 'Color', fillcol)

    slr_scenarios{ii} = SLR;
    t_scenarios{ii} = tt;


end

% artifical examples
%
% slr_rise = [0.8,2.5,3.6]; %mm slr by 2100
% slr_rise_high = 4.5* ((t-t(1))/100) - 1*((t-t(1))/100).^2;
% %slr_rise_mid  =  2.5* ((t-t(1))/100) - 0.4*((t-t(1))/100).^2;
% slr_rise_mid  = 1.2*((t-t(1))/100)+   2.8*((t-t(1))/100).^2;
% slr_rise_low = 0.8*(t-t(1))/100;
%
% plot(ax(1), slr_rise_high, t, 'k--', 'linewidth', 1.5)
% plot(ax(1), slr_rise_mid, t, 'k--', 'linewidth', 1.5)
% plot(ax(1), slr_rise_low, t, 'k--', 'linewidth', 1.5)
% slr_scenarios = {slr_rise_low;slr_rise_mid;slr_rise_high};
% t_scenarios = {t',t',t'};

%
% compute enhancement along the slr contours
%
for is = 1:3

    % get curve info
    slr_scenario = cell2mat(slr_scenarios(is));
    t_scenario   = cell2mat(t_scenarios(is));


    %init storage
    AER_scenario = nan(1,length(t_scenario));
    AER_upper_scenario = nan(1,length(t_scenario));
    AER_lower_scenario = nan(1,length(t_scenario));

    % for each time, get the enhancement
    for it = 1:length(t_scenario)

        %find the nearest time in t, which indexes AER
        [~,tidx] = min(abs(t - t_scenario(it)));

        %find hte nearest xx point to current SLR
        [~,xxidx] = min(abs(xx - slr_scenario(it)));


        AER_scenario(it) = AER(tidx, xxidx);
        AER_upper_scenario(it) = AER_upper(tidx,xxidx);
        AER_lower_scenario(it) = AER_lower(tidx,xxidx);

        %         issig_pt05_scenario(it) = issig_pt05(it,idx);
        %         issig_pt01_scenario(it) = issig_pt01(it,idx);

    end

    %adjust for zero division

    AER_upper_scenario(isinf(AER_upper_scenario)) = 1e6;
    AER_upper_scenario(isnan(AER_upper_scenario)) = 1e6;
    AER_lower_scenario(AER_lower_scenario==0) = 1e-6;
    AER_upper_scenario(AER_upper_scenario==0) = 1e-6;
    AER_lower_scenario(isinf(AER_lower_scenario)) = 1e-6; % divided by zero error
    AER_lower_scenario(isnan(AER_lower_scenario)) = 1e-6; % zero divided by zero error




    % fill the uncertainty
    idxkeep = slr_scenario > 0.1;
    xf = [t_scenario(idxkeep);flip(t_scenario(idxkeep))];

    yf = [log10(smooth(AER_lower_scenario(idxkeep))); flip(log10(smooth(AER_upper_scenario(idxkeep)))) ];
    %     plot(ax(is+1), t, log10(smooth(AER_upper_scenario)), 'r--', 'linewidth', 1.5)
    %     plot(ax(is+1), t, log10(smooth(AER_lower_scenario)), 'r--', 'linewidth', 1.5)

    fill(ax(is+1), xf,yf, fillcol, 'FaceAlpha',0.2, 'LineStyle','none');
    plot(ax(is+1), t_scenario(idxkeep), log10(smooth(AER_scenario(idxkeep))), 'k', 'linewidth', 1.5, 'Color', fillcol)


    plot(ax(is+1), t_scenario, zeros(size(t_scenario)), 'k', 'linewidth' ,1.2, 'Color', 0.5*[1,1,1])

    %    drawnow; pause

end


% tidy
ax(2).YLim = [-1,1];
ax(3).YLim = [-1,1];
ax(4).YLim = [-1,1];
ax(4).XLabel.String = 'time (years)';
for is = 1:3
    ax(is+1).XLim = [10,100];
    ax(is+1).YLabel.String = 'AER';
    ax(is+1).YLabel.Position(1) =2;
    ax(is+1).XTick = 20:20:100;
    ax(is+1).YTick = -1:1;
    ax(is+1).YTickLabel = {'10^{-1}', '10^0', '10^1'};

end
for is = 1:2
    ax(is+1).XTick = ax(4).XTick;
    ax(is+1).XTickLabel = {};

end

