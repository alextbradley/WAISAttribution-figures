% Make figure 7 of the ms, showing the difference between the calibrated
% and un-calibrated data
%
% Note that the script to make figure 6 (make_figure6.m) must be run prior
% to running this script, as the output from the former is used herein.
%
% 27/02/23 ATB (aleey@bas.ac.uk), MIT licence
%
%
% Preliminaries
%
addpath('plottools')
gendata = 1; %set to 1 to pass thru the gendata loop
fs = 13; %plot fontsize


%
% load in wavi and mitgcm data
%
if gendata
    ss_wavi = load('data/WAVI-ensemble-data.mat');
    ss_wavi = ss_wavi.ss;
    ss_mit  = load('data/MITgcm-ensemble-data.mat');
    ss_mit  = ss_mit.ss;
end

%
% run info
%
gammas      = 1:5; %1:5 correspond to 8:12 * 1e-3
gammas_act  = 8:12; %what do these gamma value actually mean
ensembles   = 1:2; %1: anthro trend, 2: no trend
members     = 1:20;
tshow       = 0:100; %time values to show SLR at
timeslices  = [0,25,50,75,100]; %calibration times

%length of arrays for conveniences
lg = length(gammas);
le = length(ensembles);
lm = length(members);
lt = length(tshow);
ltc = length(timeslices);

%
% Constants
%
rhoi = 918;  %ice density
rhow = 1028; %water density
dx   = 1000;
dy   = 1000; %grid resolution

%
% get the bathymetry
%
% compute the ice bed
fpath = strcat('data/ATTR_00000/outfile.nc');
bed   = ncread(fpath, 'b', [1, 1, 1], [Inf, Inf, 1]);
float_thick = abs(rhow/rhoi *bed); %thickness at which floatation occurs

%
% Loop thru an get SLR for each ensemble, each member, for each gamma,
% at each time, as well as calibration data
%
slr_data = struct;
Dbar     = ones(le,lm,lg); %set the calibration co-efficient to unity 

for ie = 1:le
    for im = 1:lm
        for ig = 1:lg

            % compute the vaf and slr as a function of time
            hh = ss_wavi(ig,ie,im).h; %ice thickness
            idx = hh > float_thick;
            dh = hh - float_thick;
            dh(~idx) = 0;               %mask anything below float thick to 0
            vaf = squeeze(sum(sum(dh,2),1)*dx*dy); %vaf as a function of time
            slr = (vaf(1) - vaf)/ 395 / 1e9; %SLR in mm

            for it = 1:lt
                %get the index this t point corresponds to
                [~,tidx] = min(abs(ss_wavi(ig,ie,im).t - tshow(it)));

                slr_data(ie,im,ig,it).slr = slr(tidx);
            end %end loop over time show points
 

        end %end loop over gamma values
        
    end %end loop over members
end %end loop over ensembles

%% get slr curve for each ensemble and member
%choose constnats
sigma_m = 10;
sigma_g = 1;
mu = 10;
x = linspace(-1,4,1e3);

pslr_data = struct;
count = 1;

clf;
for ie = 1:le 
    for im =1:lm
        %subplot(le,10, count); hold on; box on
        for it = 1:lt
            
            slrs = [slr_data(ie,im,:,it).slr];
            D = squeeze(Dbar(ie,im,:));
            pslr = get_pslr(x, slrs, gammas_act, D', sigma_m, sigma_g, mu);
            pslr_data(ie,im,it).pslr = pslr;
             
           %  plot(x,pslr ); 
           %  ylim([0,3])
            
        end
        count = count + 1 %, drawnow; pause 
    end

end

%% work out the means of distributions
mean_pdfs = nan(le,lt,length(x));
dx = diff(x); dx = dx(1);
for it = 10:lt
for ie = 1:le
    ens_mean = nan(size(x));
    for ix = 1:length(x)
        vals = nan(1,lm); %store all values at this particular x
        for im = 1:lm
            vals(im) = pslr_data(ie,im,it).pslr(ix);
        end

        ens_mean(ix) = median(vals);
        kde = fitdist(vals','kernel');
        ens_mean(ix) = mean(kde);
    end
    ens_mean = ens_mean / sum(ens_mean) /dx;
    mean_pdfs(ie,it, :) = ens_mean;
end
it

end

%% load in the figure 6 data
fig6data = load('figure6-out.mat');


%% make the first panel
figure(1); clf;
dp = squeeze(mean_pdfs(1,11:end,:) - mean_pdfs(2,11:end,:));
t = tshow(11:end); 
dp(abs(dp)<1e-4) = 0;
p = imagesc(x, t, dp);
clim([-.3,.301]);
set(gca, 'YDir', 'normal');
colormap(cmocean('balance'));
c = colorbar;
c.Label.String = 'anthropogenic enhancement';
c.Label.FontSize = fs+2;
ax = gca;
ax.XLim = [-0.3, 3];
ax.YLim = [10,100];
ax.XTick = 0:3;
ax.FontSize = fs;
ax.FontName = 'GillSans';
ax.YLabel.String = 'time (years)';
ax.XLabel.String = 'sea level rise (mm)';
ax.XLabel.FontSize = fs+2;
ax.YLabel.FontSize = fs+2;

fig = gcf; fig.Position(3:4) = [560,420];

%% make the second panel
% pdf_anth = squeeze(mean_pdfs(1,11:end,:));
% pdf_nat = squeeze(mean_pdfs(2,11:end,:));
% ratio = (pdf_anth - pdf_nat) ./ pdf_nat;
% 
% ratio = pdf_anth ./ pdf_nat;
% cmap = cmocean('curl', 100);
% cmap = cmap(11:end-10,:);
% 
% % anthro enhancement first
% figure(2); clf;
% 
% ratio_pos = ratio;
% ratio_pos(ratio<0)= nan;
% 
% p = imagesc(x, t,log10(ratio_pos));
% set(p, 'AlphaData', ~isnan(ratio_pos))
% clim([-2,2]);
% set(gca, 'YDir', 'normal');
% ax(1) = gca;
% colormap(ax(1), cmap(41:end,:));
% 
% % natural enhancement
% ax(2) = axes();
% ratio_neg = ratio;
% ratio_neg(ratio > 0)= nan;
% p = imagesc(x, t,log10(abs(ratio_neg)));
% %clim(ax(2),[-2,2]);
% set(p, 'AlphaData', ~isnan(ratio_neg))
% colormap(ax(2), (cmap(1:40,:)));
% 
% 
% ax(2).Visible = 'off';
% ax(2).XLim = ax(1).XLim;
% ax(2).YLim = ax(1).YLim;
% 
% 
% % make the colorbar;
% axc = axes();
% c = colorbar(axc);
% axc.Colormap = cmap;
% axc.Visible = 'off';
% 
% % tidy up
% for i = 1:2
%     ax(i).XLim = [-0.3, 3];
%     ax(i).YLim = [10,100];
% end


%% make the second panel
pdf_anth = squeeze(mean_pdfs(1,11:end,:));
pdf_nat = squeeze(mean_pdfs(2,11:end,:));
ratio = (pdf_anth - pdf_nat) ./ pdf_nat;
%ratio(pdf_nat < 1e-5) = nan;
%ratio = pdf_anth./pdf_nat;

idx = (pdf_nat ~= 0 & pdf_anth == 0);
ratio(idx) = -100;
%ratio(isinf(ratio))
cmap = cmocean('curl', 100);
cmap = cmap(11:end-10,:);

figure(2); clf;

%ratio(ratio > 0) = nan;
p = imagesc(x, t,(ratio)*100);
set(p, 'AlphaData', ~isnan(ratio))
%clim([-2,2]);
set(gca, 'YDir', 'normal');
ax(1) = gca;
colormap(ax(1), cmap);
c = colorbar;
clim([-200,200])
c.Label.String = 'percentage anthropogenic enhancement';
c.Label.FontSize = fs+2;

ax = gca;
ax.XLim = [-0.3, 3];
ax.YLim = [10,100];
ax.XTick = 0:3;
ax.FontSize = fs;
ax.FontName = 'GillSans';
ax.YLabel.String = 'time (years)';
ax.XLabel.String = 'sea level rise (mm)';
ax.XLabel.FontSize = fs+2;
ax.YLabel.FontSize = fs+2;

fig = gcf; fig.Position(3:4) = [560,420];