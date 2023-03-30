
poolobj = gcp('nocreate');
if ~isempty(poolobj);  delete(poolobj); end
num_cpu=24;
poolobj = parpool('local',num_cpu);


nboot = 1000;

anth_ci_upper = nan(length(t), length(x));
anth_ci_lower = nan(length(t), length(x));
nat_ci_upper  = nan(length(t), length(x));
nat_ci_lower  = nan(length(t), length(x));
for it = 1:length(t)
    tic
    for ix = 1:length(xx)

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
