% Make supplementary figure "B" of the manuscript showing contour plots of
% melt rate at timeslices for different M values.

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
nx = 300;
ny = 50;     %grid size

%
% Run info
%
gammas        = 1:5; %1:5 correspond to M = 0.5:1.5
ensemble      = 1; %1: anthro trend, 2: no trend
member        = 5; % e.g.. Note must be same as fig 5
timeslices    = [0,25,50,75,100];
lt =  length(timeslices);
lg = length(gammas);
%% Generate data

if gendata
    ss = load('data/WAVI-ensemble-data.mat');
    s = ss.ss;
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


%% Make the plot
% 
figure(3); clf;

count = 1;
for p = 1:length(timeslices)

    %do the wavi melt rates
    for q = 1:length(gammas)
        ax2(q,2*p-1) = subplot(2*length(timeslices),length(gammas), count);
        m = squeeze(melt_wavi(q,p,:,:));
        m(m==0) = nan;
        pl = imagesc(m');
        set(pl, 'AlphaData', ~isnan(m'));
        count = count + 1;
        
        %drawnow; pause
    end

    %then do the mitgcm melt rates

    for q = 1:lg
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

fig = gcf;
fig.Position(3:4) = [1000, 500];

%
% add a colorbar
%
c = colorbar;
c.Location = 'southoutside';
c.Position(3) = 0.3;
c.FontSize = 14;
c.Label.String = 'melt rate (m/yr)';
c.FontName = 'Arial';
c.Position(2) = 0.08;
c.Position(1) = 0.5;
c.Position(4) = 0.015;


%% Make a plot of the forcing (vertical)
fig2 = figure(2); clf; hold on; box on

% load in the data
fa = load('data/forcing_anomalies.mat');
pc = fa.ss(ensemble,member).pc;
tt  = fa.ss(ensemble,member).t;

%plot initial pos and trend
plot(-500*[1,1],[0,100] ,'k--', 'linewidth', 1.2)
plot(pc(1:20:end), tt(1:20:end), 'linewidth', 1.4, 'color', [0, 33, 153]/255);

%tidy 
axs = gca;
axs.YLim = [0,100];
axs.YTick = 0:20:100;
axs.YTickLabel = {};
axs.FontSize = 14;
axs.FontName ='Arial';
axs.XLim = [-800, -200];
axs.XTick = [-800,-500,-200];
axs.XLabel.String = '$\mathcal{P}_c(t)$'; 
axs.XLabel.Interpreter = 'latex';

% add timeslices
for it = 1:lt
   % plot([-800,200], timeslices(it)*[1,1], 'k--', 'linewidth', 1.2, 'Color',0.5*[1,1,1])
end
axs.YTick = timeslices;

fig2.Position(4) = fig.Position(4); %make same height at main panel
fig2.Position(3) = 200;