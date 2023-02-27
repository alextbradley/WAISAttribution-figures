function D = get_D(mm,mw,hh)
% return the calibration coefficient as a function of mm (mitgcm melt), mw
% (wavi melt) and hh (ice thickness)

base = -918/1028 *hh;
idx = (base < -500) & (mm ~= 0) & (mw ~=0); %points counted in the comparison
dM = abs(mm - mw);
D = mean(mean(dM(idx)));