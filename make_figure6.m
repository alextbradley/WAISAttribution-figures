% Make figure 5 of the manuscript showing:
% (1) SLR as a function of time for single ensemble member for each different value of alpha.
% (2) SLR as a function of alpha at times t1, t2, ...
% (3) Likehood as a function of gamma_T as a heatmap
% (4) P(mitgcm|mu, M) 
% (5) P(M|mu)
% (6) l(M) (normalizaed product of 4 and 5)
% (g) P(SLR = x) as a function of x
%
% 24/02/23, ATB (aleey@bas.ac.uk), MIT license
% 

%
% Preliminaries 
%
addpath('plottools')
gendata = 1; %set to 1 to pass thru the gendata loop
fs = 15; %plot fontsize
figsize = [500,300]; %set plot size


%
% Constants
%
rhoi = 918;  %ice density
rhow = 1028; %water density
dx   = 1000;
dy   = 1000; %grid resolution

%
% Run info
%
gammas        = 1:5; %1:5 correspond to 8:12 * 1e-3
Ms_act        = 0.5:0.25:1.5; %what do these gamma value actually mean in terms of M
ensemble      = 1; %1: anthro trend, 2: no trend
member        = 5; %
sz            = [5,1]; %size of ensemble
ss = load('data/WAVI-ensemble-data.mat');
s = ss.ss;


%
% Initialize storage
%
if gendata
    vaf  = cell(sz);
    t    = cell(sz); %times corresponding to vaf

    % compute the ice bed
    fpath = strcat('data/ATTR_00000/outfile.nc');
    bed   = ncread(fpath, 'b', [1, 1, 1], [Inf, Inf, 1]); %bed topo
    float_thick = abs(rhow/rhoi *bed);

    for i = 1:sz(1) %alpha value
        for j = 1:sz(2) %member number
            gamma_idx = gammas(i);
            member_idx = member(j);

            % get VAF
            hh = s(gamma_idx,ensemble,member_idx).h; %ice thickness
            idx = hh > float_thick;
            dh = hh - float_thick;
            dh(~idx) = 0;
            vv = sum(sum(dh,2),1)*dx*dy;
            vaf{i,j} = squeeze(vv);
            tt = s(gamma_idx,ensemble,member_idx).t;
            t{i,j} = tt;

        end
    end
end %end gendata flag

Mcts = linspace(min(Ms_act), max(Ms_act), 2e2); %cts form of M


%% Make plot 1: SLR as a function of time
tshow = [20,40,60,80, 100]; %times to show in b
colmapb = cmocean('haline',length(tshow)+1); %colourmap in b
colmap = flipud(cmocean('matter',sz(1) + 1)); %colourmap in a
positions = [   0.07 , 0.17, 0.4, 0.5; 
                0.55, 0.17, 0.4, 0.78];


fig(1) = figure(1); clf;
fig(1).Position(3:4) = figsize;
hold on; box on

% add tshow lines first
for it = 1:length(tshow)
    plot(tshow(it)*[1,1], [-5,5], '--', 'linewidth', 1.5, 'color', colmapb(it,:))
end

% add slr curves
for i = 1:sz(1)
    for j = 1:sz(2)

        tt = cell2mat(t(i,j));
        vv = cell2mat(vaf(i,j));
        volume_change = vv(1) - vv;
        SLR = volume_change / 395 / 1e9; %SLR in mm
        p = plot(tt, smooth(SLR,5 ), 'linewidth', 2, 'Color', colmap(i,:));

    end
end
ax1 = gca;
ax1.XLim = [0,100];
ax1.XLabel.String = 'time (years)';
ax1.FontSize = fs;
ax1.YLabel.String = '\Delta SLR (mm)';
ax1.YLim = [-1,4];
ax1.FontName = 'GillSans';

c = colorbar;
c.Colormap = (colmap(1:end-1,:));
c.Ticks = 0.1:0.2:0.9;
c.TickLabels = {'0.5', '0.75','1', '1.25','1.5'};
c.Label.String = '$M$';
c.Label.Interpreter = 'latex';
c.FontSize = fs;
c.Label.FontSize = fs;


%% Make plot 2: SLR as a function of M
fig(2) = figure(2); clf; hold on; box on;
fig(2).Position(3:4) = figsize;
%colmapb = parula(length(tshow)+1);

for it = 1:length(tshow)
    [~,idxt] = min(abs(tt - tshow(it))); %index of entries at 100 yrs (should be final entry)
    SLR_at_this_time = zeros(1,sz(1));
    for ig = 1:sz(1) %loop over gamma values
        vv = cell2mat(vaf(ig,1));
        volume_change = vv(1) - vv;
        SLR_at_this_time(ig) = volume_change(idxt) / 395 / 1e9; %SLR in mm

    end %end loop over gamma values
    SLR_times{it} = SLR_at_this_time;

    plot(Ms_act, (SLR_at_this_time),'color', colmapb(it,:), 'linewidth', 2)

%     % test the get_gamma_x script
%     x = 0.2; %text value
%     gamma_x = get_gamma_x(x, SLR_at_this_time, gammas_act);
%     plot([min(gammas_act), gamma_x], x*[1,1], '--', 'linewidth', 2)
%     1
%     drawnow ; pause

end %end loop over time pts
ax2 = gca;
ax2.XLabel.String = '$M$';
ax2.XLabel.Interpreter = 'latex';
ax2.YLabel.String = '\Delta SLR (mm)';
ax2.YLim = [-.5,4];
ax2.YTick = 0:4;
ax2.XTick = Ms_act;
ax2.FontSize = fs; 
ax2.FontName ='GillSans';
ax2.YTick = -1:4;


% %
% % create plot of forcing
% %
% fa = load('data/forcing_anomalies.mat');
% pc = fa.ss(ensembles,members).pc;
% tt  = fa.ss(ensembles,members).t;
% axs = axes(); hold on; box on
% plot([0,100], -500*[1,1], 'k--', 'linewidth', 1.2)
% plot(tt(1:20:end), pc(1:20:end), 'linewidth', 1.4, 'color', [0, 33, 153]/255);
% 
% axs.Position = [ax(1).Position(1),ax(1).Position(2) + ax(1).Position(4) + 0.01, ax(1).Position(3), ax(2).Position(4) - ax(1).Position(4)];
% axs.XLim = [0,100];
% axs.XTick = 0:20:100;
% axs.XTickLabel = {};
% axs.FontSize = 14;
% axs.FontName ='GillSans';
% axs.YLim = [-800, -200];
% axs.YTick = [-800,-500,-200];
% axs.YLabel.String = '$P(t;\mathcal{F})$'; 
% axs.YLabel.Interpreter = 'latex';

%% Generate data for the likelihood 

nx = 300;
ny = 50;
timeslices = [0,25,50,75,100];
gamma_idx = 1:5; %1,2,3,4,5 correspond to 8,9,10,11,12 *1e-3
lt = length(timeslices);
lg = length(gamma_idx);

if gendata
    ss_mit = load('data/MITgcm-ensemble-data.mat');
    s_mit =  ss_mit.ss;
    melt_wavi = nan(lg,lt,nx,ny);
    melt_mit  = nan(lg,lt,nx,ny);
    D = nan(lg,lt);

    for it = 1:lt %timeslice value
        for ig = 1:lg %gamma value 
            % get the wavi melt rates
            tt = s(ig,ensemble,member).t;
            s(ig,ensemble,member).fid;
            [~,idx] = min(abs(tt - timeslices(it))); %index of time entry matching
            
            mw = s(ig,ensemble,member).m(:,:,idx); %ice model melt rate
            melt_wavi(ig,it,:,:) = mw;

            %get the mitgcm melt rates
            mm = s_mit(ig,ensemble,member,it).m;
            melt_mit(ig,it,:,:) = mm;

            %get the calibration coefficient
            hh = s(ig,ensemble,member).h(:,:,idx); %ice thickness at this point
            D(ig,it) = get_D(mm,mw,hh);


        end
    end
end %end gendata flag


%% Make plot 3: Dbar values
fig(3) = figure(3); clf; fig(3).Position(3:4) = figsize;
Dt = D'; %take transpose so that times go down
DD = [Dt; mean(Dt)];
h = heatmap(DD); 
%surf(DD); view([0,90])
Dbar = mean(Dt); %taking means in time (first axis of D')

h.CellLabelFormat = '%.0f';
ax3 =gca;
ax3.YDisplayLabels = {'0', '25', '50', '75', '100', ''};
ax3.XDisplayLabels = {'0.5', '0.75', '1', '1.25', '1.5'};
cl = ax3.ColorLimits;
ax3.ColorLimits = [0, 50];
cmap =  flipud(cmocean('ice',100));
cmap = cmap(1:end-20,:);
ax3.Colormap = cmap;
ax3.FontName = 'GillSans';
ax3.FontSize = fs;
ax3.XLabel = '$M$';
ax3.YLabel = 'time (years)';

%put the x axis at top (ignore warning)
axp = struct(ax3);       %you will get a warning
axp.Axes.XAxisLocation = 'top';
axp.Axes.XLabel.Interpreter = 'latex';
%% Make plot 4: exp(-Dbar^2 / 2sigma_m^2)
% plot of likelihood contribution from Dbar
fig(4) = figure(4); clf; hold on; box on;
fig(4).Position(3:4) = figsize;

sigma_m = [10]; %option to do more if you want
dm = diff(Ms_act); dm = dm(1); %grid spacing in M
cmap = flipud(cmocean('amp', length(sigma_m) + 2));
cmap = zeros(10,3);

for ism = 1:length(sigma_m)
    yy = 1/sqrt(2*pi*sigma_m(ism)^2) * exp (-Dbar .^2 /2/sigma_m(ism)^2);
    yy = yy/(sum(yy)*dm); %normalize (mass lost outside range)

   
    %do in a continuous way
    f = fit(Ms_act',yy','SmoothingSpline','SmoothingParam', 1);
    density_calibrate = f(Mcts);% ff(ff > 1) = 1;
    plot(Mcts, density_calibrate, 'linewidth', 2, 'Color',cmap(ism+1,:))
        % add the discrete points
    plot(Ms_act, yy,'ko', 'linewidth', 2, 'Color',cmap(ism+1,:),'markerfacecolor', 'w', 'markersize', 7)

    %plot(gammas_act, yy,'o', 'markeredgecolor', 'k', 'markerfaceColor',cmap(ism+1,:))
    

end

fig = gcf;
fig.Position(3:4) = [500,300];

ax4 = gca; 
ax4.FontSize = fs;
ax4.XTick = Ms_act;
ax4.YLim = [0,2.5];
ax4.YTick = 0:2;
ax4.XLim = [0.5, 1.5];
ax4.FontName = 'GillSans';
ax4.YLabel.String = '$P(\dot{m}_{\mathrm{ocean}}|M)$';
ax4.XLabel.String = '$M$';
ax4.YLabel.Interpreter = 'latex';
ax4.XLabel.Interpreter = 'latex';

%%  Make plot 5: exp(-(M-mu)^2 / 2sigma_mu^2)
fig(5) = figure(5); clf; hold on; box on;
mu = 1;
sigma_mu = [0.1]; %1, 2.5, 5];
cmap = flipud(cmocean('tempo', length(sigma_mu)+2));
cmap = zeros(10,3); %make black
dx = diff(Mcts); dx = dx(1);
for ig = 1:length(sigma_mu)
    density_prior = 1/sqrt(2*pi*sigma_mu(ig)^2)*exp(-(Mcts - mu).^2 /2/sigma_mu(ig)^2);
    density_prior = density_prior/(sum(density_prior) *dx);
    % add discrete points

    plot(Mcts, density_prior, 'linewidth', 2, 'Color',cmap(ig+1,:));

    %add discrete pts
    plot(Ms_act,  1/sqrt(2*pi*sigma_mu(ig)^2)*exp(-(Ms_act - mu).^2 /2/sigma_mu(ig)^2), 'ko', 'linewidth', 2,'markerfacecolor','w');
    
end

ax5 = gca; 
ax5.FontSize = fs;
ax5.XTick = Ms_act;
ax5.FontName = 'GillSans';
ax5.YLim = [0,4.2];
ax5.FontName = 'GillSans';
ax5.YLabel.String = '$P(M|\mu)$';
ax5.XLabel.String = '$M$';
ax5.YLabel.Interpreter = 'latex';
ax5.XLabel.Interpreter = 'latex';
fig(5).Position(3:4) = figsize;
%% Make plot 6: l(M)
sigma_mu = 0.1; %single values now
sigma_m = 10;  


Ldisc = 1/sqrt(2*pi*sigma_mu^2)*exp(-(Ms_act - mu).^2 /2/sigma_mu^2) .*  ...
   ( 1/sqrt(2*pi*sigma_m^2) * exp (-Dbar .^2 /2/sigma_m^2));
Ldisc = Ldisc / (sum(Ldisc)*dm); %normalize

%or do the cts version 
L = (density_calibrate') .* density_prior;
L = L/(sum(L)*dx); %normalize

fig(6) = figure(6);clf; hold on; box on; fig(6).Position(3:4) = figsize;

%add the prior first
plot(Mcts, density_prior, 'k--', 'linewidth', 2)
plot(Mcts, L, 'k', 'linewidth', 2);
%plot(Ms_act, Ldisc, 'ko', 'linewidth', 2, 'MarkerFaceColor', 'w', 'MarkerSize', 8);

ax6 = gca; 
ax6.FontSize = fs;
ax6.XTick = Ms_act;
ax6.FontName = 'GillSans';

ax6.FontName = 'GillSans';
ax6.YLabel.String = '$P(M|\mu)$';
ax6.XLabel.String = '$M$';
ax6.YLabel.Interpreter = 'latex';
ax6.XLabel.Interpreter = 'latex';

%% Make panel 7: sea level rise predictions
fig(7) = figure(7); clf; hold on; box on
fig(7).Position(3:4) = figsize;
slr = linspace(0,4,1e2); %target slr values
ds = diff(slr);
ds = ds(1);

for it = 1:length(tshow)
    pslr = get_pslr(slr, cell2mat(SLR_times(it)), Ms_act, Dbar, sigma_m, sigma_mu, mu);
    pslr = pslr / (sum(pslr)*ds);

    plot(slr, pslr, 'linewidth', 2, 'Color', colmapb(it,:));
end 



ax7 = gca; 
ax7.FontSize = fs;
ax7.XTick = 0:3;
ax7.XLim = [0,3.5];
ax7.FontName = 'GillSans';
ax7.YLim = [0,3];
ax7.YTick = 0:3;

ax7.YLabel.String = '\Delta SLR (mm)';
ax7.XLabel.String = '$P(\Delta SLR| \mathcal{F}_i)$';
ax7.XLabel.Interpreter = 'latex';

