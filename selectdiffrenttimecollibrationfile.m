function selectdiffrenttimecollibrationfile()
cd '/home/ttflinac/user/mflorian/PolariX/FL2/time_calibration'
temp_file='/home/ttflinac/user/mflorian/PolariX/FL2/time_calibration/time_calib_9FL2XTDS.mat';
filename=uigetfile;

display(filename)
load(filename)
save(temp_file, 'amplitude_XTDS', 'timecal_fspixel', 'timecal_fspixel_err',...
    'timecal_streak', 'time_res', 'timestamp')

cd /System/Volumes/Data/home/ttflinac/user/mirian/application-for-Polarix-main






end