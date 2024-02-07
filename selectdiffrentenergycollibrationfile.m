function selectdiffrentenergycollibrationfile()
%cd '/home/ttflinac/user/mflorian/PolariX/FL2/time_calibration'
cd '/home/ttflinac/user/mflorian/PolariX/FL2/energy_calibration';
filename=uigetfile;

display(filename)
load(filename)
save('/home/ttflinac/user/mflorian/PolariX/FL2/energy_calibration/erg_calib_9FL2XTDS.mat', 'ergcal', 'ergcal_err', 'timestamp')
cd /System/Volumes/Data/home/ttflinac/user/mirian/application-for-Polarix-main


end