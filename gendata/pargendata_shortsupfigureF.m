% Generate data for supplmentary figure F of the manuscript, showing the
% AER as a function of SLR and time for different values of sigma_g (called
% sigma_P in the ms) and sigma_m (called sigma_L in the ms). Here, we
% return a nm x ng structure named data_out (where nm is the number of
% sigma_m values and ng the number of sigma_g values), each of which
% contains fields:
% x: slr values
% t: time values
% mean_pdfs: (size 2 x nx x nt) pdfs from the distributions. First 'row'
% for anthro and second 'row' for counterfactual

%% Preliminaries
%
% Parallel info
%
poolobj = gcp('nocreate');
if ~isempty(poolobj);  delete(poolobj); end
num_cpu=24;
num_cpu=2;
poolobj = parpool('local',num_cpu);

%
% Preliminaries
%
addpath('..')
addpath('../plottools')
gendata = 1; %set to 1 to pass thru the gendata loop

%
% Bayesian parameters
%
sigma_ms = [1,5,10,20];
sigma_gs = [0.01, 0.1,0.2, 0.5];
mu = 1.25;
data_out = struct();

%
% load in wavi and mitgcm data
%
if gendata
    ss_wavi = load('../data/WAVI-ensemble-data.mat');
    ss_wavi = ss_wavi.ss;
    ss_mit  = load('../data/MITgcm-ensemble-data.mat');
    ss_mit  = ss_mit.ss;
end

%
% run info
%
Ms          = 1:5; %indices of M
Ms_act      = 0.5:0.25:1.5; %what do these gamma value actually mean
ensembles   = 1:2; %1: anthro trend, 2: no trend
members     = 1:40;
tshow       = 10:100; %time values to show SLR at
timeslices  = [0,25,50,75,100]; %calibration times


%length of arrays for conveniences
lg = length(Ms);
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
x = linspace(-1,4,1e4); %output points for SLR (i.e. target values). Need very fine resolution to capture situations with little sea level rise as a fucntion of M

%
% get the bathymetry
%
% compute the ice bed
fpath = strcat('../data/ATTR_00000_outfile.nc');
bed   = ncread(fpath, 'b', [1, 1, 1], [Inf, Inf, 1]);
float_thick = abs(rhow/rhoi *bed); %thickness at which floatation occurs

%% Get the SLR data
%
% Loop thru an get SLR for each ensemble, each member, for each gamma,
% at each time, as well as calibration data
%
slr_data = struct;
Dbar     = nan(le,lm,lg); %for storing the mean calibration coefficients
allD     = nan(le,lm,lg,ltc); %for storing all calibration coefficients
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

%% Loop over different sigma_P and sigma_L values
for isigma_m = 1:length(sigma_ms)
    for isigma_g = 1:length(sigma_gs)
        pslr_data = struct;
        count = 1;

        for ie = 1:le
            for im =1:lm
                %subplot(le,10, count); hold on; box on
                parfor it = 1:lt

                    slrs = [slr_data(ie,im,:,it).slr];
                    D = squeeze(Dbar(ie,im,:));
                    pslr = get_pslr(x, slrs, Ms_act, D', sigma_ms(isigma_m), sigma_gs(isigma_g), mu);
                    pslr_data(ie,im,it).pslr = pslr;


                end
                count = count + 1 %, drawnow; pause
            end

        end




        % store all values of slr
        vals = nan(lt,le,lm,length(x)); %for each x, store all associated values
        for it = 1:lt
            for ie = 1:2
                for ix = 1:length(x)
                    for im = 1:lm
                        vals(it,ie,im,ix) = pslr_data(ie,im,it).pslr(ix);
                    end
                end
            end
            it

        end
        mean_pdfs = squeeze(mean(vals,3));

        %store the pdfs
        data_out(isigma_m, isigma_g).mean_pdfs = mean_pdfs;
        data_out(isigma_m, isigma_g).x = x;
        data_out(isigma_m, isigma_g).t = tshow;


    end %end loop over sigma_g
end %end loop over sigma_m



%% save the output
t = tshow;
save('shortshortsupfigure-Fdata.mat',"data_out");
