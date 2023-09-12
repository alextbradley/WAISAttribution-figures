% Make supplementary figure F, showing the AER as a function of SLR and
% time for different values of sigma_L and sigma_P. We make each of the
% panels individually and assemble them externally.

addpath('plottools');

%
% parameter values (should match info in structure!)
%
sigma_Ls = [5,10,20];
sigma_Ps = [0.1,0.2, 0.5];

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
    50
    100]; %intervals of colormap (choose middle number to match clims in a)

cmap = interp1(x/100,T,linspace(0,1,255));

for i = 2:4
    for j = 2:4

        fig = figure(1); clf; ax(1) = gca;
        fig.Position(3:4) = [560, 420];
        % subsample and smoothing parameters
        nxs = 10; %how finely to subsample x
        nsmooth = nxs*3;
        t       = data_out(i,j).t;
        x       = data_out(i,j).x;

        %extract the distributions
        anth_mean = data_out(i,j).mean_pdfs(:,1,:);
        nat_mean = data_out(i,j).mean_pdfs(:,2,:);

        %subsample stuff
        anth_mean_subsamp = anth_mean(:,1:nxs:end);
        nat_mean_subsamp  = nat_mean(:,1:nxs:end);
        xx = x(1:nxs:end);

        %smooth stuff
        anth_mean_subsamp_smooth = nan(size(anth_mean_subsamp));
        nat_mean_subsamp_smooth = nan(size(nat_mean_subsamp));

        %loop over time points
        for it = 1:length(t)
            anth_mean_subsamp_smooth(it,:) = smooth(squeeze(anth_mean_subsamp(it,:)), nsmooth)';
            nat_mean_subsamp_smooth(it,:) = smooth(squeeze(nat_mean_subsamp(it,:)), nsmooth)';
        end

        AER_smooth = nan(length(t), length(xx));
        for it = 1:length(t)

            % smooth distributions
            %store
            AER_smooth(it,:) = anth_mean_subsamp_smooth(it,:)./nat_mean_subsamp_smooth(it,:);

        end

        %AER(isinf(AER)) = 1e4;
        %AER_smooth = smooth2a(AER_smooth, 2,5);
        p = imagesc(ax(1), xx, t, log10(AER_smooth));

        set(p, 'AlphaData', ~isnan(AER_smooth));
        set(gca, 'YDir', 'normal');
        %c = colorbar(ax(1));
        %c.Ticks = -1:0.5:1.5;
        %c.TickLabels = {'10^{-1}','10^{-0.5}' '10^{0}', '10^{0.5}','10^{1}', '>10^{1.5}'};
        colormap(ax(1), cmap)
        ax(1).CLim = [-1,1];
        ax(1).XLim = [-0.3, 4];
        ax(1).YLim = [10,100];
        %ax(1).YLabel.String = 'time (years)';
        %ax(1).XLabel.String = 'SLR (mm)';
        ax(1).YTick = 20:20:100;
        ax(1).YTickLabel = {};
        ax(1).XTick = 0:4;
        ax(1).XTickLabel = {};
        title(['\sigma_L = '  num2str(sigma_Ls(i)), ', \sigma_P = ' num2str(sigma_Ps(j))])

        % save it
        exportgraphics(gcf, strcat('figures-short/raw/supfigF/panel_',num2str(i),'_',num2str(j),'.pdf'), 'ContentType', 'Vector');
    end
end

%% make plots of the likelihood functions
% (1): 1/sqrt(2*pi*sigma_L^2)*exp(-D^2 / 2 sigma_L^2) as a function of D
fig2 = figure(2); clf; hold on; box on; 
fig2.Position(3:4) = [560, 320];
sigma_Ls = [5,10,20];
D = linspace(0,20);
colmap = parula(6);
clear legendinfo
for i = 1:3
   plot(D,1/sqrt(2*pi*sigma_Ls(i)^2)*exp(-D.^2 / 2 /sigma_Ls(i)^2), 'linewidth', 1.5, 'color', colmap(i,:) );
    legendinfo{i} = ['$\sigma_L = ' num2str(sigma_Ls(i)), '$ (m/yr)'];
   
end
axs = gca;
legend(legendinfo, 'Interpreter','latex', 'FontSize', 14);
axs.FontName = 'Arial';
axs.FontSize = 14;
axs.XLabel.String = 'D (m/yr)';
axs.YLabel.String = '$(1/\sqrt{2\pi\sigma_L^2})\exp[-D^2 /(2\sigma_L^2)]$';
axs.YLabel.Interpreter = 'latex';

%%
%(2) alpha / sqrt(2*pi*sigma_P^2) * exp(-(M - mu)^2 / 2/sigma_P^2) as a
%function of M, with mu = 1.25

fig3 = figure(3); clf; hold on; box on; 
fig3.Position(3:4) = [560, 320];
sigma_Ps = [ 0.1,0.2, 0.5];
mu = 1.25;
M = linspace(0.5,1.5, 200);
dM = diff(M); dM = dM(1);
colmap = cool(5);
clear legendinfo
for j = 1:3
    pp = 1 / sqrt(2*pi*sigma_Ps(j)^2) * exp(-(M - mu).^2 / 2/sigma_Ps(j)^2);
    pp = pp/(sum(pp)*dM); %normalize the distn
    plot(M,pp, 'linewidth', 1.5, 'color', colmap(j+1,:) );
    legendinfo{j} = ['$\sigma_P = ' num2str(sigma_Ps(j)), '$'];

end
axs = gca;
legend(legendinfo, 'Interpreter','latex', 'FontSize', 14, 'location', 'northwest');
axs.FontName = 'Arial';
axs.FontSize = 14;
axs.XLabel.String = '$M$';
axs.YLabel.String = '$(\alpha/\sqrt{2\pi\sigma_P^2})\exp[-(M - \mu)^2 /(2\sigma_P^2)]$';
axs.YLabel.Interpreter = 'latex';axs.XLabel.Interpreter = 'latex';