
poolobj = gcp('nocreate');
if ~isempty(poolobj);  delete(poolobj); end
num_cpu=24;
poolobj = parpool('local',num_cpu);

pdata = load('shortfigure-3data.mat');

nxs = 10; %how finely to subsample x
nsmooth = nxs*3;
t            = pdata.t;

%extract the distributions 
anth_all           = pdata.vals(:,1,:,:);
nat_all            = pdata.vals(:,2,:,:);
anth_mean          = squeeze(mean(anth_all,3)); %take mean over ensemble members
nat_mean           = squeeze(mean(nat_all,3)); %take mean over ensemble members

%subsample stuff
anth_all_subsamp  = squeeze(anth_all(:,:,:,1:nxs:end));
nat_all_subsamp   = squeeze(nat_all(:,:,:,1:nxs:end));
anth_mean_subsamp = anth_mean(:,1:nxs:end);
nat_mean_subsamp  = nat_mean(:,1:nxs:end);

%smooth stuff
anth_all_subsamp_smooth  = nan(size(anth_all_subsamp));
nat_all_subsamp_smooth   = nan(size(nat_all_subsamp));
anth_mean_subsamp_smooth = nan(size(anth_mean_subsamp));
nat_mean_subsamp_smooth = nan(size(nat_mean_subsamp));

for it = 1:length(t)
    for im = 1:40
        anth_all_subsamp_smooth(it,im,:) = smooth(squeeze(anth_all_subsamp(it,im,:)), nsmooth)';
    end

end

xx                 = pdata.x(1:nxs:end);


nboot = 1000;

anth_ci_upper = nan(length(t), length(xx));
anth_ci_lower = nan(length(t), length(xx));
nat_ci_upper  = nan(length(t), length(xx));
nat_ci_lower  = nan(length(t), length(xx));
for it = 1:length(t)
    tic
    parfor ix = 1:length(xx)

        %anthro
        vals = squeeze(anth_all_subsamp(it,:,ix));
        ci = bootci(nboot,@mean,vals);
        anth_ci_upper(it,ix) = ci(1);
        anth_ci_lower(it,ix) = ci(2);

        %counterfactual
        vals = squeeze(nat_all_subsamp(it,:,ix));
        ci = bootci(nboot,@mean,vals);
        nat_ci_upper(it,ix) = ci(1);
        nat_ci_lower(it,ix) = ci(2);
        
        if mod(ix,100) == 0
            ix
        end

    end
    toc
    it
end



figure(1); clf; 
p = imagesc(xx, t, log10(ci_upper));
