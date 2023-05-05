% Make figure 9 of the ms, showing the anthropogenic enhancement plot for
% different sigma_m and sigma_g.

% 27/02/23, ATB (aleey@bas.ac.uk), MIT licence


%
% Preliminaries
%
addpath('plottools')
gendata = 1; %set to 1 to pass thru the gendata loop
fs = 13; %plot fontsize

%
%choose constnats
%
sigma_gs = [0.5,1, 2.5, 5];
sigma_ms = [5,10,15,20];
mu = 10;
x = linspace(-1,4,1e3);
fig8_data = struct; %initialize a structure to store the pslr and ensemble means as a function of sigma_g and sigma_m

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

%% for each sigma_m, sigma_g, get slr curve for each ensemble and members
for isg = 1:length(sigma_gs)
    for ism = 1:length(sigma_ms)
        sigma_g = sigma_gs(isg);
        sigma_m = sigma_ms(ism);

        pslr_data = struct;
        count = 1;

        %
        % get the pslr data
        %
        for ie = 1:le
            for im =1:lm
                for it = 1:lt

                    slrs = [slr_data(ie,im,:,it).slr];
                    D = squeeze(Dbar(ie,im,:));
                    pslr = get_pslr(x, slrs, gammas_act, D', sigma_m, sigma_g, mu);
                    pslr_data(ie,im,it).pslr = pslr;

                end
                count = count + 1 %
            end

        end

        %
        % work out the means of distributions
        %
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
        %store this data
        fig8_data(isg,ism).mean_pdfs = mean_pdfs;
        fig8_data(isg,ism).pslr = pslr_data;
    end %end loop over sigma_m
end %end loop over sigma_g

%% Make the plot
figure(1); clf; 
count = 1;
for isg = 1:length(sigma_gs)
    for ism = 1:length(sigma_ms)
        ax(isg,ism) = subplot(length(sigma_gs), length(sigma_ms), count);
        mean_pdfs = fig8_data(isg, ism).mean_pdfs;
        t = 10:100;
        anthro_enhance = squeeze(mean_pdfs(1,11:end,:) - mean_pdfs(2,11:end,:));
        p = imagesc(x, t, smooth2a(anthro_enhance,10,10));
        clim([-.3,.3]);
        set(gca, 'YDir', 'normal');
        colormap(cmocean('balance'));
     %  title(sprintf('s_g = %.2f, s_m = %.2f', sigma_gs(isg), sigma_ms(ism)))
        
        if isg == length(sigma_gs)
            ax(isg,ism).XTick = 0:3;
            ax(isg,ism).XLabel.String = 'sea level rise (mm)';
        else
            ax(isg,ism).XTick = [];
        end
        
        if ism == 1
            ax(isg,ism).YTick = 20:20:100;
            ax(isg,ism).YLabel.String = 'time (years)';

        else
            ax(isg,ism).YTick = [];
        end

        ax(isg,ism).XLim = [-0.3, 3];
        ax(isg,ism).FontSize = 13;
        ax(isg,ism).FontName = 'GillSans';

        count = count + 1;
    end
end
fig = gcf;
fig.Position(3:4) = [1040,730];
