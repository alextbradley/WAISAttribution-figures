% Make figure 9 of the ms, showing the difference between the calibrated
% and un-calibrated data
%
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
% load in data for calibrated (figure 8) and uncalibrated (figure 9)
%
fig8_calibrated = load("figure8-data.mat");
fig9_uncalibrated = load("figure9-data.mat");


% 
% sub-sample variables
% 
calib_pdf = squeeze(fig8_calibrated.mean_pdfs(:,1, 1:50:end));
uncalib_pdf = squeeze(fig9_uncalibrated.mean_pdfs(:,1, 1:50:end));
xx           = fig8_calibrated.x(1:50:end);
t            = fig8_calibrated.t;

%
% compute change
%
dp_calib = calib_pdf - uncalib_pdf;
dp_calib = smooth2a(dp_calib, 2,3);

% 
% make plot
%
p = imagesc(xx, t, dp_calib);
set(gca, 'YDir', 'normal');
colormap(cmocean('delta'))
clim([-.3,.301]);

c = colorbar;
c.Label.String = 'uncalibrated - calibrated';
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



