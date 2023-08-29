% Make figure 3 of the manuscript showing (a) calibrated distributions of
% sea level rise through time for both anthropogenic and counterfactual
% ensembles and (b)--(d) summary statistics of the two distributions with
% and without calibration. 
%
% 27/02/23, ATB (aleey@bas.ac.uk), MIT licence
%

%% Preliminaries
% load in the data
pdata = load('data/shortfigure-3data.mat');

addpath('plottools');
colmap = nan(2,3);
colmap(1,:) = [255,152,51]/255;  %anthro
colmap(2,:) = [0,153, 153]/255; %counter

% for poster
% colmap(2,:) = [7,54,80]/255; %dark blue
% colmap(1,:) = [248,200,44]/255; %yellow

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
fig.Position(3:4) = [1080, 400];


%% Make (a): joy plots of distributions
% choose output times
nxs = 10; %how finely to subsample x
nsmooth = nxs*3;
t            = pdata.t;
anth_mean_pdf_allt = squeeze(pdata.mean_pdfs(:,1, 1:nxs:end));
nat_mean_pdf_allt = squeeze(pdata.mean_pdfs(:,2, 1:nxs:end));

tplot_idx    = 1:10:length(t); %indices to show on joy plot

% sub sample the mean pdfs and is_significant matrices
anth_mean_pdf = squeeze(pdata.mean_pdfs(tplot_idx,1, 1:nxs:end));
nat_mean_pdf = squeeze(pdata.mean_pdfs(tplot_idx,2, 1:nxs:end));
xx           = pdata.x(1:nxs:end);
dx = diff(xx); dx = dx(1);
t_plot       = t(tplot_idx);

for i = flip(1:length(t_plot))

    baseline = t_plot(i);

    % plot the nat pdf
    height = 18;
    xf = [xx, flip(xx)];

    % smooth, renormalize and scale distributions
    ysnat = smooth(nat_mean_pdf(i,:), nsmooth)'; 
    ysnat = ysnat / sum(ysnat * dx);

    ysanth = smooth(anth_mean_pdf(i,:), nsmooth)';
    ysanth = ysanth / sum(ysanth * dx);
    
    % compute medians
    runsum = cumsum(ysanth *dx); [~,medidx] = min(abs(runsum - 0.5));
    mean_anth = xx(medidx);
    runsum = cumsum(ysnat *dx); [~,medidx] = min(abs(runsum - 0.5));
    mean_nat = xx(medidx);

    ymax = max([ysnat,ysanth]);
    ysnatsc = ysnat*height/ymax;
    ysanthsc = ysanth*height/ymax;


    %fill natural
    yf = [baseline*ones(size(xx)), flip(baseline*ones(size(xx))+ysnatsc)];
    fill(ax(1), xf, yf,colmap(2,:), 'FaceAlpha', 0.5);
    plot(ax(1), xx, baseline*ones(size(xx))+ysnatsc, 'k', 'linewidth', 1, 'color', 0.2*[1,1,1])

    %fill anthro
    yf = [baseline*ones(size(xx)), flip(baseline*ones(size(xx))+ysanthsc)];

    fill(ax(1), xf, yf, colmap(1,:), 'FaceAlpha', 0.5);
    plot(ax(1), xx, baseline*ones(size(xx))+ysanthsc, 'k', 'linewidth', 1, 'color', 0.2*[1,1,1])

    % add the mean
    plot(ax(1),mean_nat, baseline , 'o', 'MarkerFaceColor',colmap(2,:), 'MarkerEdgeColor','k')
    plot(ax(1),mean_anth, baseline , 'o', 'MarkerFaceColor',colmap(1,:), 'MarkerEdgeColor','k', 'MarkerSize',7.5)

end

ax(1).XLim = [-0.3,4];
%ax(1).View = [18,50];

ax(1).YTick = t_plot;
ax(1).YLim = [min(t_plot), max(t_plot)+height+1];


ax(1).XLabel.String = 'sea level rise (mm)';
ax(1).YLabel.String = 'time (years)';
box(ax(1), 'off');
ax(1).Color = 'none';
shg

%% Make (b): showing ensemble mean and std
hold on; box on;

skewnesses = nan(length(t), 2);
kurtosises = nan(length(t), 2);
means      = nan(length(t), 2);


lt = length(t);

for it = 1:lt

    anth_pdf = smooth(squeeze(anth_mean_pdf_allt(it,:)), nsmooth);anth_pdf = anth_pdf / sum(anth_pdf*dx);
    nat_pdf  = smooth(squeeze(nat_mean_pdf_allt(it,:)), nsmooth); nat_pdf = nat_pdf / sum(nat_pdf*dx);

    %compute means
    means(it,2)      = sum(nat_pdf' .* xx) *dx;
    means(it,1)      = sum(anth_pdf' .* xx) *dx;

    %compute medians
    runsum = cumsum(anth_pdf' *dx); [~,medidx] = min(abs(runsum - 0.5));
    medians(it,1) = xx(medidx);

    runsum = cumsum(nat_pdf' *dx); [~,medidx] = min(abs(runsum - 0.5));
    medians(it,2) = xx(medidx);

    %compute fourth moments
    fourth_moments(it,1) = sum(anth_pdf' .* (xx - means(it,1)).^4 * dx);
    fourth_moments(it,2) = sum(nat_pdf' .* (xx - means(it,2)).^4 * dx);

    %compute variance
    variance(it,1) = sum(anth_pdf' .* (xx - means(it,1)).^2 * dx);
    variance(it,2) = sum(nat_pdf' .* (xx - means(it,2)).^2 * dx);

    %compute kurtosis
    kurtosises(it,1) = fourth_moments(it,1)/ variance(it,1)^2;
    kurtosises(it,2) = fourth_moments(it,2)/ variance(it,2)^2;

    %compute skewnesses
    skewnesses(it,1) = sum(anth_pdf' .* ((xx - means(it,1))/sqrt(variance(it,1))).^3 * dx);
    skewnesses(it,2) = sum(nat_pdf' .* ((xx - means(it,2))/sqrt(variance(it,2))).^3 * dx);

end


% plot means
plot(ax(2), t, smooth(medians(:,1),5), 'linewidth', 2, 'color', colmap(1,:));
plot(ax(2), t, smooth(medians(:,2),5), 'linewidth', 2, 'color', colmap(2,:));

% plot skewness in fourth plot
plot(ax(3), t, smooth(skewnesses(:,1),5), 'linewidth', 2, 'color', colmap(1,:));
plot(ax(3), t, smooth(skewnesses(:,2),5), 'linewidth', 2, 'color', colmap(2,:));


% plot kurtosis in third plot
plot(ax(4), t, smooth(kurtosises(:,1),5), 'linewidth', 2, 'color', colmap(1,:));
plot(ax(4), t, smooth(kurtosises(:,2),5), 'linewidth', 2, 'color', colmap(2,:));


%ax(3).YLim = [0,6];

%% Repeat for the uncalibrated data 
pdata_uncalib = load('data/shortfigure3-uncalibdata.mat');
t_uncalib            = pdata_uncalib.t;
anth_mean_pdf_allt_uncalib = squeeze(pdata_uncalib.mean_pdfs(:,1, 1:nxs:end));
nat_mean_pdf_allt_uncalib = squeeze(pdata_uncalib.mean_pdfs(:,2, 1:nxs:end));
xx_uncalib          = pdata_uncalib.x(1:nxs:end);

for it = 1:lt

    anth_pdf = smooth(squeeze(anth_mean_pdf_allt_uncalib(it,:)), nsmooth);
    anth_pdf = anth_pdf / sum(anth_pdf*dx);
    nat_pdf  = smooth(squeeze(nat_mean_pdf_allt_uncalib(it,:)), nsmooth);
    nat_pdf = nat_pdf / sum(nat_pdf*dx);

    %compute means
    means(it,2)      = sum(nat_pdf' .* xx) *dx;
    means(it,1)      = sum(anth_pdf' .* xx) *dx;

    %compute medians
    runsum = cumsum(anth_pdf' *dx); [~,medidx] = min(abs(runsum - 0.5));
    medians(it,1) = xx(medidx);

    runsum = cumsum(nat_pdf' *dx); [~,medidx] = min(abs(runsum - 0.5));
    medians(it,2) = xx(medidx);

    %compute fourth moments
    fourth_moments(it,1) = sum(anth_pdf' .* (xx - means(it,1)).^4 * dx);
    fourth_moments(it,2) = sum(nat_pdf' .* (xx - means(it,2)).^4 * dx);

    %compute variance
    variance(it,1) = sum(anth_pdf' .* (xx - means(it,1)).^2 * dx);
    variance(it,2) = sum(nat_pdf' .* (xx - means(it,2)).^2 * dx);

    %compute kurtosis
    kurtosises(it,1) = fourth_moments(it,1)/ variance(it,1)^2;
    kurtosises(it,2) = fourth_moments(it,2)/ variance(it,2)^2;

    %compute skewnesses
    skewnesses(it,1) = sum(anth_pdf' .* ((xx - means(it,1))/sqrt(variance(it,1))).^3 * dx);
    skewnesses(it,2) = sum(nat_pdf' .* ((xx - means(it,2))/sqrt(variance(it,2))).^3 * dx);

end


% plot means
plot(ax(2), t, smooth(medians(:,1),5),'--',  'linewidth', 2, 'color', colmap(1,:));
plot(ax(2), t, smooth(medians(:,2),5), '--', 'linewidth', 2, 'color', colmap(2,:));

% plot skewness in fourth plot
plot(ax(3), t, smooth(skewnesses(:,1),5), '--','linewidth', 2, 'color', colmap(1,:));
plot(ax(3), t, smooth(skewnesses(:,2),5), '--','linewidth', 2, 'color', colmap(2,:));


% plot kurtosis in third plot
plot(ax(4), t, smooth(kurtosises(:,1),5),'--', 'linewidth', 2, 'color', colmap(1,:));
plot(ax(4), t, smooth(kurtosises(:,2),5), '--','linewidth', 2, 'color', colmap(2,:));

%% tidy
for i = 2:4
    ax(i).XLim = [10,100];
    ax(i).YLabel.Position(1) = 3.5;
end

for i = 2:3
    ax(i).XTick = ax(4).XTick;
    ax(i).XTickLabel = {};
end

ax(4).XLabel.String = 'time (years)';
ax(4).YLabel.String = 'kurtosis';
ax(3).YLabel.String = 'skewness';
ax(2).YLabel.String = 'median (mm)';
    