

function get_time_duration(firstphase, secondphase)
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
        %charge(jj)=doocsread(addr_charge);

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

        pause(0.11)

    end


    img_mean = squeeze(mean(img_sig, 3));

    %%%

    % ccharge_7FL2XTDS=,

    a=1
   [ length_y, length_x]=size(img_mean );
    %y_axis, y_axis= collibrataxis(length_y, length_x)
    imag=hlc_clean_image( squeeze(img_mean(a:end,:)));

    % time
    tmp_profile         = mean(imag);
    [x_com(i), x_var(i), x_fwhm(i), x_axis, x_profile] = get_profile_stats(tmp_profile, timecal_fspixel);


end

i=2;
tem=doocswrite(addr_onbeam,0)
display(' phase was set')



for jj=1:nshot
    ddd_read            = doocsread([addr_cam, 'IMAGE_EXT_ZMQ']);
    img_sig(:,:,jj)     = ddd_read.data.val_val;

    ts(jj)=ddd_read.timestamp;
 %   charge(jj)=doocsread(addr_charge);

    % compared it with last measurement
    if jj > 1
        while ts(jj) == ts(jj-1)
            ddd_read            = doocsread([addr_cam, 'IMAGE_EXT_ZMQ']);
            ts(jj)=ddd_read.timestamp;
            img_sig(:,:,jj)       = ddd_read.data.val_val;
            display( [' - (', num2str(jj), ') same data ... wait ...']);
            pause(0.1)
        end
    end

    pause(0.1)

end
tem=doocswrite(addr_onbeam,1)

img_mean = squeeze(mean(img_sig, 3));

%%%

% ccharge_7FL2XTDS=,


[length_y, length_x]=size(img_sig);
%y_axis, y_axis= collibrataxis(length_y, length_x)
imag=hlc_clean_image( squeeze(img_sig(a:end,:)));

% time
 tmp_profile         = mean(imag);
    [x_com(i), x_var(i), x_fwhm(i), x_axis, x_profile]  = get_profile_stats(tmp_profile, timecal_fspixel);

%%
volt=[-1,0, 1]

figure(122321)
plot(volt,x_fwhm, '-.o')
xlim([-1.2 1.2])
xticks([-1, 0, 1])
xticklabels({'-v','0','v'})

