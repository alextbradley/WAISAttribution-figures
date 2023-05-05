% Generate data to produce figure 10. Produces variables:
% (1) mean_pdfs, size (nt x ne x nx), the mean probability distribution in
% the uncalibrated case

% nt: number of time output points
% ne: number of ensemlbes (here 2)
% nm: number of ensemble members (here 20)
% nx: number of tagret slr points 
%
% NB: this script runs in a parallel fashion.
% 
% Files save these, alongside time output points, in a file figure8-out.mat
%
% NB this script takes ~O(10s) of minutes to run.
%
% ATB (aleey@bas.ac.uk), 17/3/23. MIT licence.

% 
% Parallel info
%
poolobj = gcp('nocreate');
if ~isempty(poolobj);  delete(poolobj); end
num_cpu=24;
poolobj = parpool('local',num_cpu);

%
% load in wavi and mitgcm data
%
gendata = 1;
if gendata
    ss_wavi = load('data/WAVI-ensemble-data.mat');
    ss_wavi = ss_wavi.ss;
    ss_mit  = load('data/MITgcm-ensemble-data.mat');
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
sigma_gs = [0.5,1, 2.5, 5];
sigma_ms = [5,10,15,20];
mu = 1;
x = linspace(-1,4,1e4); %output points for SLR (i.e. target values). Need very fine resolution to capture situations with little sea level rise as a fucntion of M 

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
Dbar     = ones(le,lm,lg); %for storing the mean calibration coefficients
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

%% get slr curve for each realization of forcing
varsig_data = struct;
count = 1;
for isg = 1:length(sigma_gs)
    for ism = 1:length(sigma_ms)
        sigma_g = sigma_gs(isg);
        sigma_m = sigma_ms(ism);
        pslr_data = struct;
        count = 1;
        vals = nan(lt,le,lm,length(x)); %for each x, store all associated values


        for ie = 1:le
            for im =1:lm
                %subplot(le,10, count); hold on; box on
                parfor it = 1:lt

                    slrs = [slr_data(ie,im,:,it).slr];
                    D = squeeze(Dbar(ie,im,:));
                    pslr = get_pslr(x, slrs, Ms_act, D', sigma_m, sigma_g, mu);
                    pslr_data(ie,im,it).pslr = pslr;

                    vals(it,ie,im,:) = pslr;

                end
                count = count + 1; %, drawnow; pause
            end

        end
        mean_pdfs = squeeze(mean(vals,3));
        varsig_data(isg, ism).mean_pdfs = mean_pdfs;
	count = count + 1
    end
end



%% save the output for use in figure 10
t = tshow;
save('figure10-data.mat', "varsig_data", "t", "x");
