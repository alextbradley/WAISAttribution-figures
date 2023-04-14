% Make supplementary figure A, showing the proportion of simulations that
% have retreated as a function of time.
%
% 31/03/23, ATB (aleey@bas.ac.uk)


%% Preliminaries
% load in the data
%data = load('data/WAVI-ensemble-data.mat');
wavdat = data.ss;
gendata = 1; %pass through loop to generate
%% Figure setup
fig = figure(1); clf;
fig.Color = 'w';
fig.Position(3:4) = [560, 340];
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Arial';
hold(ax, 'on')
box(ax, 'on')


% colours
colmap = nan(2,3);
colmap(1,:) = [255,152,51]/255;  %anthro
colmap(2,:) = [0,153, 153]/255; %counter

%% Process the data
if gendata
    % compute the slr through time for each simulation
    SLRs = nan(5,2,40,101);
    for iM = 1:5
        for ie  = 1:2
            for im = 1:40
                hh = wavdat(iM,ie,im).h; %ice thickness
                idx = hh > float_thick;
                dh = hh - float_thick;
                dh(~idx) = 0;
                vv = sum(sum(dh,2),1)*dx*dy;
                vv = squeeze(vv);
                volume_change = vv(1) - vv;
                SLR = volume_change / 395 / 1e9; %SLR in mm
                SLRs(iM,ie,im,:) = SLR;
                tt = wavdat(iM,ie,im).t;   %same in each case!

            end
        end
    end

    % turn into a retreat proportion
    propret = nan(2,5,101);
    threshret = 0.2;
    for ie = 1:2
        for iM = 1:5
            for it = 1:101

                SLR_members = squeeze(SLRs(iM,ie,:,it));
                propret(ie,iM,it) = sum(SLR_members > threshret)/length(SLR_members);

            end
        end
    end

    proprete = squeeze(mean(propret,2));
end %end gendata flag
%% Plot stuff
for i =[2,1]
    plot(tt, proprete(i,:), 'Color', colmap(i,:), 'LineWidth',1.75)
end

% add regression?
propanth = proprete(1,:); idx = find(propanth>0,1, 'First');idx = idx-1;
propanth = propanth(idx:end); tanth = tt(idx:end);
p = polyfit(tanth,propanth,1); 
plot(tanth, p(1)*tanth+p(2), '--', 'Color', colmap(1,:), 'LineWidth', 1.5);

propnat = proprete(2,:); idx = find(propnat>0,1, 'First');idx = idx-1;
propnat = propnat(idx:end); tanth = tt(idx:end);
pp = polyfit(tanth,propnat,1); 
plot(tanth, pp(1)*tanth+pp(2), '--', 'Color', colmap(2,:), 'LineWidth', 1.5);


ax.XLabel.String = 'time (years)';
ax.XLabel.Position(2) = -0.08;
ax.YLabel.String = 'proportion retreated';
ax.YLabel.Position(1)  = -10;
legend({'Counterfactual', 'Anthropogenic'}, 'location', 'Northwest');
ax.YLim = [0,1];
ax.Position(2) = 0.13;
shg

