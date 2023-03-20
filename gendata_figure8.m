% Generate data to produce figure 8. Produces variables:
% (1) mean_pdfs, size (nt x ne x nx), the mean probability distribution 
% (2) vals, size (nt x ne x nm x nx), value of probabilities.
% (3) is_significant_ptX (nt x ne x nx) boolean matrices, showing whether
% the set of probabilities are significantly different from one another at
% significance level X.
% nt: number of time output points
% ne: number of ensemlbes (here 2)
% nm: number of ensemble members (here 20)
% nx: number of tagret slr points 
% 
% Files save these, alongside time output points, in a file figure8-out.mat
%
% NB this script takes ~O(10s) of minutes to run.
%
% ATB (aleey@bas.ac.uk), 17/3/23. MIT licence.


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
sigma_m = 10;
sigma_g = 0.1;
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

%% get slr curve for each realization of forcing

pslr_data = struct;
count = 1;

clf;
for ie = 1:le 
    for im =1:lm
        %subplot(le,10, count); hold on; box on
        for it = 1:lt
            
            slrs = [slr_data(ie,im,:,it).slr];
            D = squeeze(Dbar(ie,im,:));
            pslr = get_pslr(x, slrs, Ms_act, D', sigma_m, sigma_g, mu);
            pslr_data(ie,im,it).pslr = pslr;
             
           %  plot(x,pslr ); 
           %  ylim([0,3])
            
        end
        count = count + 1 %, drawnow; pause 
    end

end




%% store all values of slr
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

%% get the significance curves
is_significant_pt1 = zeros(lt,length(x));
is_significant_pt05 = zeros(lt,length(x));
is_significant_pt01 = zeros(lt,length(x)); %store significance at different levels
for it = 1:lt
    for ix = 1:length(x)
    
        v1 = squeeze(vals(it,1,:,ix)); %all anthro members at this time and x
        v2 = squeeze(vals(it,2,:,ix)); %all non anthro members
        [~,h] = ranksum(v1, v2, 'Alpha', 0.1);
        is_significant_pt1(it,ix) =  h;
        [~,h] = ranksum(v1, v2, 'Alpha', 0.05);
        is_significant_pt05(it,ix) =  h;
        [~,h] = ranksum(v1, v2, 'Alpha', 0.01);
        is_significant_pt01(it,ix) =  h;
    end
    it
end

%% save the output for use in figure 7
t = tshow;
save('figure8-out.mat', "mean_pdfs", "vals", "t", "is_significant_pt1", "is_significant_pt05","is_significant_pt01", "x");
