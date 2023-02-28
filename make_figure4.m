% Make figure 4 of the manuscript showing:
% (a) SLR as a function of time for single ensemble member for each different value of alpha.
% (b) SLR as a function of alpha at times t1, t2, ...
% (c) Table of contour plots of melt rate at timeslices for different alpha values. Likelihood values and mean likelihood labelled.
% (d) Likehood as a function of gamma_T as a heatmap
% (e) P(mitgcm|mu, gamma) (for different values of sigma_m?)
% (f) P(gamma|mu) (for different values of sigma_a)
% (g) P(SLR = x) as a function of x
% (nb this script produces (c) last)
%
% 24/02/23, ATB (aleey@bas.ac.uk), MIT license
% 

%
% Preliminaries
%
addpath('plottools')
gendata = 1; %set to 1 to pass thru the gendata loop
fs = 13; %plot fontsize

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
gammas_act    = 8:12; %what do these gamma value actually mean
ensembles     = 1; %1: anthro trend, 2: no trend
members       = 5; %5 is good exampel
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
            member_idx = members(j);

            % get VAF
            hh = s(gamma_idx,ensembles,member_idx).h; %ice thickness
            idx = hh > float_thick;
            dh = hh - float_thick;
            dh(~idx) = 0;
            vv = sum(sum(dh,2),1)*dx*dy;
            vaf{i,j} = squeeze(vv);
            tt = s(gamma_idx,ensembles,member_idx).t;
            t{i,j} = tt;

        end
    end
end %end gendata flag

%% make panels a and b
% Plot data
tshow = [20,40,60,80, 100]; %times to show in b
colmapb = cmocean('haline',length(tshow)+1); %colourmap in b
colmap = flipud(cmocean('matter',sz(1) + 1)); %colourmap in a
positions = [   0.07 , 0.17, 0.4, 0.5; 
                0.55, 0.17, 0.4, 0.78];

%
% Make (a)
%
fig1 = figure(1); clf; fig = gcf; fig.Position(3:4) = [1000,250];
ax(1) = subplot('Position', positions(1,:));
hold on; box on


% add tshow lines first
for it = 1:length(tshow)
    plot(tshow(it)*[1,1], [-5,5], '--', 'linewidth', 1.2, 'color', colmapb(it,:))
end
for i = 1:sz(1)
    for j = 1:sz(2)


        tt = cell2mat(t(i,j));
        vv = cell2mat(vaf(i,j));
        volume_change = vv(1) - vv;
        SLR = volume_change / 395 / 1e9; %SLR in mm
        p = plot(tt, smooth(SLR,5 ), 'linewidth', 2, 'Color', colmap(i,:));

    end
end
ax(1).XLim = [0,100];
ax(1).XLabel.String = 'time (years)';
ax(1).FontSize = 12;
ax(1).YLabel.String = 'SLR (mm)';
ax(1).YLim = [-1,4];
c = colorbar;
c.Colormap = (colmap(1:end-1,:));
c.Ticks = 0.1:0.2:0.9;
c.TickLabels = {'8', '9','10', '11','12'};
c.Label.String = '$\gamma_T~\times 10^{3}$';
c.Label.Interpreter = 'latex';
c.FontSize = 14;
c.Label.FontSize = 16;


%
% make (b)
%

ax(2) = subplot('Position',positions(2,:)); hold on; box on
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
    plot(gammas_act, (SLR_at_this_time),'color', colmapb(it,:), 'linewidth', 2)


%     % test the get_gamma_x script
%     x = 0.2; %text value
%     gamma_x = get_gamma_x(x, SLR_at_this_time, gammas_act);
%     plot([min(gammas_act), gamma_x], x*[1,1], '--', 'linewidth', 2)
%     1
%     drawnow ; pause

end %end loop over time pts

ax(2).XLabel.String = '$\gamma_T~\times 10^{3}$';
ax(2).XLabel.Interpreter = 'latex';
ax(2).YLabel.String = 'SLR (mm)';
ax(2).YLim = [-.5,4];
ax(2).YTick = 0:4;
ax(2).XTick = 8:12;

%
% tidy stuff
%
for i =1:2
    ax(i).FontSize = 14;
    ax(i).FontName ='GillSans';
end
ax(1).Position(4) = 0.5; %make smaller to accomodate forcing plot
ax(1).YTick = -1:2:3;

%
% create plot of forcing
%
fa = load('data/forcing_anomalies.mat');
pc = fa.ss(ensembles,members).pc;
tt  = fa.ss(ensembles,members).t;
axs = axes(); hold on; box on
plot([0,100], -500*[1,1], 'k--', 'linewidth', 1.2)
plot(tt(1:20:end), pc(1:20:end), 'linewidth', 1.4, 'color', [0, 33, 153]/255);

axs.Position = [ax(1).Position(1),ax(1).Position(2) + ax(1).Position(4) + 0.01, ax(1).Position(3), ax(2).Position(4) - ax(1).Position(4)];
axs.XLim = [0,100];
axs.XTick = 0:20:100;
axs.XTickLabel = {};
axs.FontSize = 14;
axs.FontName ='GillSans';
axs.YLim = [-800, -200];
axs.YTick = [-800,-500,-200];
axs.YLabel.String = '$P(t;\mathcal{F})$'; 
axs.YLabel.Interpreter = 'latex';

%% generate data for the rest of the plots
gendata_c = 1; 

nx = 300;
ny = 50;
timeslices = [0,25,50,75,100];
gamma_idx = 1:5; %1,2,3,4,5 correspond to 8,9,10,11,12 *1e-3
lt = length(timeslices);
lg = length(gamma_idx);
ensemble = ensembles;
ensemble =1;

member = members; %same as above?
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


%% make panel d
figure(4); clf; 
Dt = D'; %take transpose so that times go down
DD = [Dt; mean(Dt)];h = heatmap(DD); %take the 
Dbar = mean(Dt); %taking means in time (first axis of D')

 h.CellLabelFormat = '%.0f';
ax =gca;
ax.YDisplayLabels = {'0', '25', '50', '75', '100', ''};
ax.XDisplayLabels = {'8', '9', '10', '11', '12'};
cl = ax.ColorLimits;
ax.ColorLimits = [0, cl(2)];
cmap =  flipud(cmocean('ice',100));
cmap = cmap(1:end-20,:);
ax.Colormap = cmap;
ax.FontName = 'GillSans';
ax.FontSize = fs;
fig = gcf;
fig.Position(3:4) = [500,300];

%% make panel e
% plot of likelihood contribution from Dbar
figure(5); clf; hold on; box on;

%create a spline fit to this 

xx = linspace(min(gammas_act), max(gammas_act), 5e3);
sigma_m = [5,10,15,20]; 

cmap = flipud(cmocean('amp', length(sigma_m) + 2));
for ism = 1:length(sigma_m)
    yy = 1/sqrt(2*pi*sigma_m(ism)^2) * exp (-Dbar .^2 /2/sigma_m(ism)^2);


    yy = yy/sum(yy); %normalize (grid spacing is 1)

    f = fit(gammas_act',yy','Gauss1');
    ff = f(xx); ff(ff > 1) = 1;
    plot(xx, ff, 'linewidth', 2, 'Color',cmap(ism+1,:))
    %plot(gammas_act, yy,'o', 'markeredgecolor', 'k', 'markerfaceColor',cmap(ism+1,:))
    
end
fig = gcf;
fig.Position(3:4) = [500,300];

ax = gca; 
ax.FontSize = fs;
ax.XTick = gammas_act;
ax.YLim = [0,1];
ax.FontName = 'GillSans';
%ax.YLabel.String = 'density';

%% make panel f
figure(6); clf; hold on; box on;
mu = 10;
sigma_g = [0.5,1, 2.5, 5];
xx = linspace(min(gammas_act), max(gammas_act));
cmap = flipud(cmocean('tempo', length(sigma_g)+2));
dx = diff(xx); dx = dx(1);
for ig = 1:length(sigma_g)
    density = 1/sqrt(2*pi*sigma_g(ig)^2)*exp(-(xx - mu).^2 /2/sigma_g(ig)^2);
    density = density/sum(density) *dx;
    plot(xx, density, 'linewidth', 2, 'Color',cmap(ig+1,:));
end
fig = gcf;
fig.Position(3:4) = [500,300];

ax = gca; 
ax.FontSize = fs;
ax.XTick = gammas_act;
ax.FontName = 'GillSans';

%% make panel g
figure(7); clf; hold on; box on
x = linspace(0,4,1e2);
dx = diff(x);
dx = dx(1);

%choose constnats
sigma_m = 10;
sigma_g = 1;
mu = 10;

for it = 1:length(tshow)

    slr = cell2mat(SLR_times(it));
    pslr = get_pslr(x, slr, gammas_act, Dbar, sigma_m, sigma_g, mu);
    pslr = pslr / (sum(pslr)*dx);

    plot(x, smooth(pslr, 3), 'linewidth', 2, 'Color', colmapb(it,:));
end 


fig = gcf;
fig.Position(3:4) = [500,300];

ax = gca; 
ax.FontSize = fs;
ax.XTick = 0:3;
ax.XLim = [0,3];
ax.FontName = 'GillSans';
ax.YLim = [0,3];
ax.YTick = 0:3;

%% make panel c 
% 
figure(3); clf;

 ncols = length(gamma_idx);
% nrows = 2*length(timeslices);
% colgap = 0.02;
% starty = 0.02;
% rowgap = 0.01;
% height = 1/(nrows+2); %height of plot
% width  = 1/(ncols+1); %width of plot
% startx = (1 -width*ncols - (ncols-1)*colgap)/2;
% positions = zeros(4, ncols, nrows);
% for p = 1:nrows
%     for q = 1:ncols
%         positions(:,q,p) = [startx + (q-1)*colgap + (q-1)*width, starty + (p-1)*height + (p-1)*rowgap, width, height];
%        
%     end
% end

count = 1;
for p = 1:length(timeslices)
    %do the wavi melt rates
    for q = 1:length(gammas)
        
        %ax2(q,p) = subplot('Position', positions(:,q,nrows - p+1));
        ax2(q,2*p-1) = subplot(2*length(timeslices),length(gammas), count);
        m = squeeze(melt_wavi(q,p,:,:));
        m(m==0) = nan;
        pl = imagesc(m');
        set(pl, 'AlphaData', ~isnan(m'));
        count = count + 1;
        
        %drawnow; pause
    end

    %then do the mitgcm melt rates

    for q = 1:ncols
        %ax2(q,p) = subplot('Position', positions(:,q,nrows - p+1));
        ax2(q,2*p) = subplot(2*length(timeslices),length(gammas), count);
        m = squeeze(melt_mit(q,p,:,:));
         m(m==0) = nan;
        pl = imagesc(m');
        set(pl, 'AlphaData', ~isnan(m'));
        count = count + 1;
               
       % drawnow;pause
    end
   
end

%
% do the tidying for them all
%
sz = size(ax2);
for i = 1:sz(1)
    for j = 1:sz(2)
        ax2(i,j).XTick = [];
        ax2(i,j).YTick = [];
        box(ax2(q,p), 'on');
        ax2(i,j).CLim = [0,100];
        ax2(i,j).XLim = [100,300];
        ax2(i,j).Colormap = cmocean('thermal');
    end
end

fig2 = gcf;
fig2.Position(3:4) = [1000, 500];

%
% add a colorbar
%
c = colorbar;
c.Position(4) = 0.4;
c.FontSize = 12;
c.Position(1) = 0.92;
c.Label.String = 'melt rate (m/yr)';
