% Make figure 8, showing the 
% (a) anthropogenic enhancement through time 
% (b) ensemble means
% (c) percentage enhacencement of anthro
% (d) percentage enhancement along some realization of sea level rise 
%
% 27/02/23, ATB (aleey@bas.ac.uk), MIT licence
% 

load('figure8-data.mat');

fig = figure(1); clf;
fig.Position(3:4) = [1200,600];

w = 0.35; 
h = 0.38;
positions = [0.09, 0.58, w, h;
             0.56, 0.58, w, h;
             0.09, 0.1, w, h;
             0.56, 0.1, w, h];

fs = 14;
% make the first panel
ax(1) = subplot('Position', positions(1,:));

% sub sample the mean pdfs and is_significant matrices
anth_mean_pdf = squeeze(mean_pdfs(:,1, 1:50:end));
nat_mean_pdf = squeeze(mean_pdfs(:,2, 1:50:end));
xx           = x(1:50:end);
issig_pt01    = is_significant_pt01(:,1:50:end);
issig_pt1    = is_significant_pt1(:,1:50:end);
issig_pt05    = is_significant_pt05(:,1:50:end);


dp = anth_mean_pdf - nat_mean_pdf;
nsx = 3; %smoothing number in x
nst = 2; %smoothing number in t
 
dp = smooth2a(dp,nst,nsx); % smoothing (second argument is x)

%t = tshow; %(11:end); 
dp(abs(dp)<1e-4) = 0;
p = imagesc(xx, t, dp);
clim([-.25,.25]);
set(gca, 'YDir', 'normal');
colormap(cmocean('balance'));
c1 = colorbar;
c1.Position(1) = 0.44;
c1.Position(3) = 0.01;
c1.Position(1) =  positions(1,1)+positions(1,3)+ 0.005;
c1.Label.String = 'anthropogenic enhancement';
c1.Label.FontSize = fs+2;
c1.Ticks = -0.2:0.1:0.2;
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
%sigcmap = (parula());
sigcmap = lines(6);
sigcmap = [ 0.9290    0.6940    0.1250;
    0, 200/255, 197/255];

nst = 3;
nsx = 12; %smoothing for these contours

contour(xx,t,smooth2a(issig_pt01, nst, nsx), [0.9,0.9],'color',sigcmap(1,:), 'linewidth', 1.5)
contour(xx,t,smooth2a(issig_pt05, nst, nsx),[0.9,0.9],'color',sigcmap(2,:), 'linewidth', 1.5)
%contour(xx,t,smooth2a(issig_pt1, nst, nsx), [0.9,0.9],'color',sigcmap(3,:), 'linewidth', 1.2)


% Make (b): anthro enhancement percentage
ax(2) = subplot('Position', positions(2,:));
 hold on; box on
ratio = anth_mean_pdf./ nat_mean_pdf;
ratio(nat_mean_pdf < 1e-7) = nan;
ratio(anth_mean_pdf < 1e-7) = nan;

cmap = cmocean('curl', 200);
cmap = cmap(68:end,:);

pdat = smooth2a(ratio, 1,1);
ratpos = ratio; ratpos(ratio<0) = nan; 
p = imagesc(xx, t,log10(pdat));
set(p, 'AlphaData', ~isnan(ratio))
colormap(ax(2), cmap)
set(gca, 'YDir', 'normal');
clim([-1,3])

% add significance contours
contour(xx,t,smooth2a(issig_pt01, nst, nsx), [0.9,0.9],'color',sigcmap(1,:), 'linewidth', 1.5)
contour(xx,t,smooth2a(issig_pt05, nst, nsx),[0.9,0.9],'color',sigcmap(2,:), 'linewidth', 1.5)

% tidy
ax(2).XLim = [-0.3, 3];
ax(2).YLim = [10,100];
ax(2).XTick = 0:3;
ax(2).FontSize = fs;
ax(2).FontName = 'GillSans';
ax(2).YLabel.String = 'time (years)';
ax(2).XLabel.String = 'sea level rise (mm)';
ax(2).XLabel.FontSize = fs+2;
ax(2).YLabel.FontSize = fs+2;

%sort out the colorbar;
c2 = colorbar;
c2.Position(1) =  positions(2,1)+positions(2,3)+ 0.005;c2.Position(3) = 0.01;
c2.Label.FontSize = fs+2;
c2.Label.Interpreter = 'latex';
c2.Label.String = '$P_{\mathrm{anth}}/P_{\mathrm{none}}$';

c2.Ticks= -1:3;
c2.TickLabels = {'10^{-1}', '10^{0}','10^{1}','10^{2}','10^{3}'};


% Make (c): ensemble mean
ax(3) = subplot('Position', positions(3,:));
hold on; box on;
colmap = nan(2,3);
colmap(2,:) = [ 0,    0.45    0.84];  %blue for no trend
colmap(1,:) = [ 0.80    0.24    0.1]; %red for anthro
dx = diff(xx); dx = dx(1);
lt = length(t);
for ie = 1:2
    std_slr = nan(1,lt);
    mean_slr  = nan(1,lt);
    for it = 1:lt
        if ie == 1
            pdf_at_this_time = squeeze(anth_mean_pdf(it,:));
        else
            pdf_at_this_time = squeeze(nat_mean_pdf(it,:));
        end
            mean_slr(it) = sum(pdf_at_this_time .* xx) *dx; %E(X) = int x f(x) dx   
            std_slr(it)  = sqrt(sum((xx - mean_slr(it)).^2 .*  pdf_at_this_time) * dx);
        
 
    end
    
    % add background with errors
    fx = [t, flip(t)]; fy = [smooth((mean_slr - std_slr),10); smooth(flip(mean_slr + std_slr), 10)];
    fill(fx, fy, colmap(ie,:), 'FaceAlpha',0.2, 'LineStyle','none')

    % add ensemble means
    plot(t, smooth(mean_slr), 'linewidth', 2, 'color', colmap(ie,:))

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


% (d) synthetic example

slr_rise = 2; %mm slr by 2100
slr = slr_rise*(t/100).^2;
plot(ax(2), slr, t, 'k', 'linewidth', 2)

% work out enhancement along this curve
ratio_line = zeros(1,lt);
for i = 1:lt
    slr_now = slr(i);
    [~,idx]= min(abs(slr_now - xx));
    ratio_line(i) = ratio(i,idx);
end


ax(4) = subplot('Position', positions(4,:));
hold on; box on;

%add background significance levels !! hard coded for slr_rise = 2
fill(ax(4), [70,100, 100, 70], [0.1,0.1,  1e3, 1e3],sigcmap(2,:), 'linestyle', 'none')
fill(ax(4), [80,100, 100, 80], [0.1,0.1,  1e3, 1e3],sigcmap(1,:), 'linestyle', 'none')

plot(t,smooth(ratio_line,9), 'k', 'linewidth',2)
set(gca, 'YScale', 'log');

ax(4).YLim = 10.^[-0.5,2.2];
ax(4).XLim = [10,100];
ax(4).FontSize = fs;
ax(4).FontName = 'GillSans';
ax(4).XLabel.String = 'time (years)';
ax(4).YLabel.Interpreter = 'latex';
ax(4).YLabel.String = '$P_{\mathrm{anth}}/P_{\mathrm{none}}$';
ax(4).XLabel.FontSize = fs+2;
ax(4).YLabel.FontSize = fs+2;

 ax(3).XTick = 20:20:100;
ax(4).XTick = 20:20:100;