

function filename=get_time_duration(firstphase, secondphase)
addr_TDSphase='FLASH.RF/LLRF.CONTROLLER/CTRL.POLARIX/SP.PHASE';
% adress of onbeam =======>
addr_onbeam='FLASH.DIAG/TIMINGINFO/FLFXTDS/ON_BEAM';
addr_cam        = 'FLASH.DIAG/CAMERA/OTR9FL2XTDS/';
load('/home/ttflinac/user/mflorian/PolariX/FL2/time_calibration/time_calib_9FL2XTDS.mat', 'timecal_fspixel', ...
    'timecal_fspixel_err', 'timecal_streak', 'time_res');
load('/home/ttflinac/user/mflorian/PolariX/FL2/energy_calibration/erg_calib_9FL2XTDS.mat', 'ergcal', 'ergcal_err')
display('data is loaded');

addpath('/System/Volumes/Data/home/ttflinac/user/mirian/FL2_funs')
a=1
nshot=10;
fontSize        = 14;
calib_x             = timecal_fspixel;
calib_y             = ergcal;
timeres_calib       = time_res;
   a=1
phase(1)= firstphase;
phase(3)=secondphase;
phase(2)=0.0
%% measuremtn and data saving
for i=1:2:3
    tem=doocswrite(addr_TDSphase,phase(i))
    display(' phase was set')



    for jj=1:nshot
        ddd_read            = doocsread([addr_cam, 'IMAGE_EXT_ZMQ']);
        img_sig(:,:,jj)     = ddd_read.data.val_val;

        ts(jj)=ddd_read.timestamp;

        % compared it with last measurement
        if jj > 1
            while ts(jj) == ts(jj-1)
                ddd_read            = doocsread([addr_cam, 'IMAGE_EXT_ZMQ']);
                img_sig(:,:,jj)       = ddd_read.data.val_val;
                ts(jj)=ddd_read.timestamp; 
                display( [' - (', num2str(jj), ') same data ... wait ...']);
                pause(0.1)
            end
        end
        
 
    imag=hlc_clean_image( squeeze(img_sig(a:end,:,jj)));

    % time
    tmp_profile         = mean(imag);
    [x_com(i, jj), x_var(i,jj), x_fwhm(i,jj), x_axis, x_profile] = get_profile_stats(tmp_profile, timecal_fspixel);
 pause(0.11)

    end
  image(i)=squeeze(mean(img_sig, 3)); % average of image 

end
%%
i=2;
tem=doocswrite(addr_onbeam,0)
display(' phase was set')



for jj=1:nshot
    ddd_read            = doocsread([addr_cam, 'IMAGE_EXT_ZMQ']);
    img_sig(:,:,jj)     = ddd_read.data.val_val;
    ts(jj)=ddd_read.timestamp;


    % compared it with last measurement
    if jj > 1
        while ts(jj) == ts(jj-1)
            ddd_read            = doocsread([addr_cam, 'IMAGE_EXT_ZMQ']);
            ts(jj)=ddd_read.timestamp;
            img_sig(:,:,jj)       = ddd_read.data.val_val;
            display( [' - (', num2str(jj), ') same data ... wait ...']);
            pause(0.1)
        end
    imag=hlc_clean_image( squeeze(img_sig(a:end,:,jj)));

    % time
    tmp_profile         = mean(imag);
    [x_com(i, jj), x_var(i,jj), x_fwhm(i,jj), x_axis, x_profile] = get_profile_stats(tmp_profile, timecal_fspixel);
 pause(0.11)

    end
 image(i)=squeeze(mean(img_sig, 3)); % average of image 
    pause(0.1)

end
%%
errotr_value=std(x_var, 2);
mean_value=mean(x_var, 2);
%%
S=[-1,0, 1]

figure(122321)
errorbar(S,mean_value,errotr_value, 'o', "MarkerSize",10,...
    "MarkerEdgeColor",'m',"MarkerFaceColor",'b', 'LineWidth',3)
xlim([-1.2 1.2])
xticks([-1, 0, 1])
xticklabels({'-S','0','S'})
xlable('S_y')
ylable('\sigma(fs)')
%%
p = polyfit(S,mean_value.^2,2)
%%
t=date;
Year=year(t);
fulder_add=['/home/ttflinac/data/PolariX/', num2str(Year), '/'];
timestamp       = datestr(clock, 'yyyy-mm-ddTHHMMSS');
save([fulder_add, timestamp, '.mat'],'timestamp' ,  'x_var', 'x_fwhm','image',  'timecal_fspixel')
output_filename=[fulder_add, timestamp, '.mat']
