% Make figure 4 of the manuscript, showing (a) the anthropogenic
% enhancement ratio and (b)--(d) example curves along this trajectory
%
% 27/02/23, ATB (aleey@bas.ac.uk), MIT licence
%

%% Preliminaries
% load in the data
%pdata = load('shortfigure3-data.mat');

addpath('plottools');

%% Plot setup
fig = figure(1); clf;
fig.Color = 'w';

ht = 0.85; %total height
gapy = 0.05;
hts = (ht - 2*gapy)/3;
starty = 0.11;

positions = [0.13, starty, 0.35, 0.85;
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

%create colormap between these
T = [colmap(2,:);
     1,1,1;
     colmap(1,:)];

x = [0
     40
     100]; %intervals of colormap (choose middle number to match clims in a)

cmap = interp1(x/100,T,linspace(0,1,255));

% colourmap for significance levels
sigcmap = [ 0.9290    0.6940    0.1250;
            0, 200/255, 197/255];
sigcmap = lines(2);
sigcmap = [231, 76, 60;
           108, 52, 131 ]/255;
%% (a) Anthropogenic enhancement and significance contours
% subsample and smoothing parameters
nxs = 10; %how finely to subsample x
nsmooth = nxs*3;
t            = pdata.t;
anth_mean_pdf_allt = squeeze(pdata.mean_pdfs(:,1, 1:nxs:end));
nat_mean_pdf_allt  = squeeze(pdata.mean_pdfs(:,2, 1:nxs:end));
anth_mean_pdf      = squeeze(pdata.mean_pdfs(tplot_idx,1, 1:nxs:end));
nat_mean_pdf       = squeeze(pdata.mean_pdfs(tplot_idx,2, 1:nxs:end));
issig_pt01         = pdata.is_significant_pt01(:,1:nxs:end);
issig_pt1          = pdata.is_significant_pt1(:,1:nxs:end);
issig_pt05         = pdata.is_significant_pt05(:,1:nxs:end);
xx                 = pdata.x(1:nxs:end);

AER = nan(length(t), length(xx));
for i = 1:length(t)

    % smooth distributions
    ysnat = smooth(nat_mean_pdf_allt(i,:), nsmooth)';
    ysant = smooth(anth_mean_pdf_allt(i,:), nsmooth)';

    %store
    AER(i,:) = ysant./ ysnat;

end

AER(isinf(AER)) = 1e4;
p = imagesc(ax(1), xx, t, log10(AER));

set(p, 'AlphaData', ~isnan(AER));
set(gca, 'YDir', 'normal');
c = colorbar(ax(1));
c.Ticks = -1:0.5:1.5;
c.TickLabels = {'10^{-1}','10^{-0.5}' '10^{0}', '10^{0.5}','10^{1}', '>10^{1.5}'};
colormap(ax(1), cmap)
ax(1).CLim = [-1,1.5];
ax(1).XLim = [-0.3, 4];
ax(1).YLim = [10,100];
ax(1).XLabel.String = 'time (years)';
ax(1).YTick = 20:20:100;

% 
% add significance contorurs
%

contour(ax(1), xx,t,smooth2a(issig_pt01, 0, nsmooth), [0.9,0.9],'color',sigcmap(1,:), 'linewidth', 1.75)
contour(ax(1), xx,t,smooth2a(issig_pt05, 0, nsmooth),[0.9,0.9],'color',sigcmap(2,:), 'linewidth', 1.75)


%% Make (b) showing different trajectories
slr_rise = [0.8,2.5,3.6]; %mm slr by 2100

slr_rise_high = 4.5* ((t-t(1))/100) - 1*((t-t(1))/100).^2;
slr_rise_mid  =  2.5* ((t-t(1))/100) - 0.4*((t-t(1))/100).^2;
slr_rise_low = 0.8*(t-t(1))/100;

plot(ax(1), slr_rise_high, t, 'k--', 'linewidth', 1.5)
plot(ax(1), slr_rise_mid, t, 'k--', 'linewidth', 1.5)
plot(ax(1), slr_rise_low, t, 'k--', 'linewidth', 1.5)

slr_scenarios = [slr_rise_low;slr_rise_mid;slr_rise_high];

for is = 1:3
    %compute the enhancement along this contour
    AER_scenario = nan(1,length(t));
    slr_scenario = slr_scenarios(is,:);
    issig_pt05_scenario = nan(1,length(t));
    issig_pt01_scenario = nan(1,length(t));

    
    for it = 1:length(t)
        [~,idx] = min(abs(xx - slr_scenario(it)));
        AER_scenario(it) = AER(it,idx);

%         issig_pt05_scenario(it) = issig_pt05(it,idx);
%         issig_pt01_scenario(it) = issig_pt01(it,idx);
        

    end
    plot(ax(is+1), t, ones(size(t)), 'k--')
    plot(ax(is+1), t, smooth(AER_scenario), 'k', 'linewidth', 1.5)
end


% tidy
ax(2).YLim = [0,2];
ax(3).YLim = [0,4];
ax(4).YLim = [0,10];
ax(4).XLabel.String = 'time (years)';
for is = 1:3
    ax(is+1).XLim = [10,100];
    ax(is+1).YLabel.String = 'AER';
    ax(is+1).YLabel.Position(1) =6;
end
for is = 1:2
    ax(is+1).XTick = ax(4).XTick; 
    ax(is+1).XTickLabel = {};
end

