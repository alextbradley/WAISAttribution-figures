% Make figure 5 of the manuscript showing:
% (a--d) pdfs of individual ensemble members at times (a) 40, (b) 60, (c)
% 80, (d) 100 years
% (e) combined pdfs

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
tshow       = [40,60,80,100]; %time values to show SLR at
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
Dbar     = nan(le,lm,lg); %for storing the calibration coefficients

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

            %
            % get the calibration values
            %
            D_here = nan(1,ltc);
            for itc = 1:ltc
                [~,tidx] = min(abs(ss_wavi(ig,ie,im).t - timeslices(itc)));

                %get the mitgcm melt rates
                m_mit = ss_mit(ig,ie,im,itc).m;
             
                %get the wavi melt rate
                m_wavi = ss_wavi(ig,ie,im).m(:,:,tidx); %ice model melt rate
                
              
                %get the calibration coefficient assoc w/ this timeslice
                hh = ss_wavi(ig,ie,im).h(:,:,tidx); %ice thickness at this point
                D_here(itc) = get_D(m_mit,m_wavi,hh); %mean over 'calibration points'
                allD(ie,im,ig,itc) = D_here(itc);
                
            end %end loop over timeslice claibration points
            Dbar(ie,im,ig) = mean((D_here)); %mean over the timeslices   

        end %end loop over gamma values
        
    end %end loop over members
end %end loop over ensembles

%% get slr curve for each
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
        count = count + 1;% drawnow; pause 
    end

end

%% work out the means of distributions
mean_pdfs = nan(le,lt,length(x));
for it = 1:lt
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
    mean_pdfs(ie,it, :) = ens_mean;
end

end

%% make (a)--(d)
colmap = lines(2); %colormap to differentiate between anthro and non
colmap(1,:) = [ 0,    0.45    0.84];
colmap(2,:) = [ 0.80    0.24    0.1];
%colmap = [0, 63, 153; 153,0,63]/255;
av = 0.12; %alpha value
dx = diff(x); dx = dx(1);
clf; hold on; box on;
for it = 1:4 %new figure for each timeslice
    subplot(2,2,it); hold on; box on;
    for ie = 1:2
        for im = 1:lm
            p = plot(x,pslr_data(ie,im,it).pslr, 'linewidth', 1.2, 'color', colmap(ie,:));
            p.Color = [colmap(ie,:), av];

        end

         %add the ensemble mean
         mpdf = squeeze(mean_pdfs(ie,it,:));
         mpdf = mpdf / sum(mpdf) /dx;
         plot(x, smooth(mpdf, 8),  'linewidth', 2, 'color', colmap(ie,:));
    end
    

    %tidy the plot
    ax(it) = gca;
    ax(it).YLim = [0,3];
    ax(it).XLim = [-0.2, 2];
    ax(it).YTick = 0:3;
    ax(it).XTick = 0:3;
    ax(it).FontName = 'GillSans';
    ax(it).FontSize = fs;
    ax(it).XLabel.String = 'sea level rise (mm)';
    ax(it).YLabel.String = 'density';
end

fig = gcf; fig.Position(3:4) = [1060,600];