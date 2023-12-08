
function erg_axis, time_axis=collibrataxis(length_y, length_x)

load('/home/ttflinac/user/mflorian/PolariX/FL2/time_calibration/time_calib_9FL2XTDS.mat', 'timecal_fspixel', ...
    'timecal_fspixel_err', 'timecal_streak', 'time_res');
load('/home/ttflinac/user/mflorian/PolariX/FL2/energy_calibration/erg_calib_9FL2XTDS.mat', 'ergcal', 'ergcal_err')
display('data is loaded');

calib_x             = timecal_fspixel;
calib_y             = ergcal;
timeres_calib       = time_res;
erg_axis            = ergcal * ((1:length_y)-length_y/2);
time_axis           = timecal_fspixel * ((1:length_x)-length_x/2);