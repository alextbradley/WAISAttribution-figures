% Make figure 8, showing the 
% (a) anthropogenic enhancement through time 
% (b) ensemble means
% (c) percentage enhacencement of anthro
% (d) percentage enhancement along some realization of sea level rise 
%
% 27/02/23, ATB (aleey@bas.ac.uk), MIT licence
% 

%load('figure8-out.mat');

fig = figure(1); clf;
fig.Position(3:4) = [1200,800];
%% make the first panel
ax(1) = subplot(2,2,1);
dp = squeeze(mean_pdfs(:,1, :) - mean_pdfs(:,2,:));
% for it = 1:lt
%     dp(it,:) = smooth(dp(it,:), 100); %apply smoothing in x
% end
ns = 100; %smoothing number
dp = smooth2a(dp,0,ns); % smoothing (second argument is x)

t = tshow; %(11:end); 
dp(abs(dp)<1e-4) = 0;
p = imagesc(x, t, dp);
clim([-.3,.3]);
set(gca, 'YDir', 'normal');
colormap(cmocean('balance'));
c = colorbar;
c.Label.String = 'anthropogenic enhancement';
c.Label.FontSize = fs+2;
ax(1).XLim = [-0.3, 3];
ax(1).YLim = [10,100];
ax(1).XTick = 0:3;
ax(1).FontSize = fs;
ax(1).FontName = 'GillSans';
ax(1).YLabel.String = 'time (years)';
ax(1).XLabel.String = 'sea level rise (mm)';
ax(1).XLabel.FontSize = fs+2;
ax(1).YLabel.FontSize = fs+2;


% 
% add contours of significance
%
hold on
sigcmap = cmocean('algae',6);
sigcmap = (parula(5));

contour(x,t,smooth2a(is_significant_pt05,1,10), [0.9,0.9],'color',sigcmap(3,:), 'linewidth', 1.2)
contour(x,t,smooth2a(is_significant_pt05,1,10), [0.9,0.9],'color',sigcmap(4,:), 'linewidth', 1.2)
contour(x,t,smooth2a(is_significant_pt01,1,10), [0.9,0.9],'color',sigcmap(5,:), 'linewidth', 1.2)


%% Make (b): anthro enhancement percentage
ax(2) = subplot(2,2,2); hold on; box on
pdf_anth = squeeze(mean_pdfs(:,1,:));
pdf_nat = squeeze(mean_pdfs(:,2,:));
ratio = (pdf_anth - pdf_nat) ./ pdf_nat;
%ratio(pdf_nat < 1e-5) = nan;
%ratio = pdf_anth./pdf_nat;

idx = (pdf_nat ~= 0 & pdf_anth == 0);
ratio(idx) = -100;
%ratio(isinf(ratio))
cmap = cmocean('curl', 100);
cmap = cmap(11:end-10,:);


p = imagesc(x, t,smooth2a(ratio,0,2)*100);
set(p, 'AlphaData', ~isnan(ratio))
%clim([-2,2]);
set(gca, 'YDir', 'normal');
colormap(ax(2), cmap);
c = colorbar;
clim([-250,250])
c.Label.String = 'percentage anthropogenic enhancement';
c.Label.FontSize = fs+2;

ax(2).XLim = [-0.3, 3];
ax(2).YLim = [10,100];
ax(2).XTick = 0:3;
ax(2).FontSize = fs;
ax(2).FontName = 'GillSans';
ax(2).YLabel.String = 'time (years)';
ax(2).XLabel.String = 'sea level rise (mm)';
ax(2).XLabel.FontSize = fs+2;
ax(2).YLabel.FontSize = fs+2;

%% Make (c): ensemble mean
ax(3) = subplot(2,2,3); hold on; box on;
colmap = nan(2,3);
colmap(2,:) = [ 0,    0.45    0.84];  %blue for no trend
colmap(1,:) = [ 0.80    0.24    0.1]; %red for anthro
dx = diff(x); dx = dx(1);
for ie = 1:2
    std_slr = nan(1,lt);
    mean_slr  = nan(1,lt);
    for it = 1:lt
        pdf_at_this_time = squeeze(mean_pdfs(it,ie, :));
        mean_slr(it) = sum(pdf_at_this_time .* x') *dx; %E(X) = int x f(x) dx   
        std_slr(it)  = sqrt(sum((x' - mean_slr(it)).^2 .*  pdf_at_this_time) * dx);
 
    end
    
    % add background with errors
    fx = [t, flip(t)]; fy = [mean_slr - std_slr, flip(mean_slr + std_slr)];
    fill(fx, fy, colmap(ie,:), 'FaceAlpha',0.2, 'LineStyle','none')

    % add ensemble means
    plot(t, mean_slr, 'linewidth', 2, 'color', colmap(ie,:))

end


ax(3).YLim = [-0.3, 2];
ax(3).XLim = [10,100];
ax(3).YTick = 0:3;
ax(3).FontSize = fs;
ax(3).FontName = 'GillSans';
ax(3).XLabel.String = 'time (years)';
ax(3).YLabel.String = 'sea level rise (mm)';
ax(3).XLabel.FontSize = fs+2;
ax(3).YLabel.FontSize = fs+2;

