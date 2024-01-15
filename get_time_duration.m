

function get_time_duration(firstphase, secondphase)
addr_TDSphase='FLASH.RF/LLRF.CONTROLLER/CTRL.POLARIX/SP.PHASE';
% adress of onbeam =======>
% addr_onbeam=

load('/home/ttflinac/user/mflorian/PolariX/FL2/time_calibration/time_calib_9FL2XTDS.mat', 'timecal_fspixel', ...
    'timecal_fspixel_err', 'timecal_streak', 'time_res');
load('/home/ttflinac/user/mflorian/PolariX/FL2/energy_calibration/erg_calib_9FL2XTDS.mat', 'ergcal', 'ergcal_err')
display('data is loaded');

addpath('/System/Volumes/Data/home/ttflinac/user/mirian/FL2_funs')

nshot=10;
fontSize        = 14;
calib_x             = timecal_fspixel;
calib_y             = ergcal;
timeres_calib       = time_res;

phase(1)= firstphase;
phase(2)=secondphase;

%% measuremtn and data saving
for i=1:2:3
    tem=doocswrite(add_phase,phase(i))
    display(' phase was set')



    for j=1:nshot
        ddd_read            = doocsread([app.addr_cam, 'IMAGE_EXT_ZMQ']);
        img(:,:,jj)     = ddd_read.data.val_val;

        ts(jj)=ddd_read.timestamp;
        charge(jj)=doocsread(addr_charge);

        % compared it with last measurement
        if jj > 1
            while ts(jj) == ts(jj-1)
                ddd_read            = doocsread([app.addr_cam, 'IMAGE_EXT_ZMQ']);
                img_sig(:,:,jj)       = ddd_read.data.val_val;
                display( [' - (', num2str(jj), ') same data ... wait ...']);
                pause(1/app.rep_rate)
            end
        end

        pause(1/app.rep_rate)

    end


    img_mean = squeeze(mean(img_sig, 3));

    %%%

    % ccharge_7FL2XTDS=,

    a=1
    length_y, length_x=size(img_sig)
    %y_axis, y_axis= collibrataxis(length_y, length_x)
    imag=hlc_clean_image( squeeze(img_sig(a:end,:)))

    % time
    tmp_profile         = img_mean;
    [x_com(i), x_var(i), x_fwhm(i), x_axis(i), x_profile(i)] = get_profile_stats(tmp_profile, timecal_fspixel);


end

i=2;
tem=doocswrite(addr_offbeam,1)
display(' phase was set')



for j=1:nshot
    ddd_read            = doocsread([app.addr_cam, 'IMAGE_EXT_ZMQ']);
    img(:,:,jj)     = ddd_read.data.val_val;

    ts(jj)=ddd_read.timestamp;
    charge(jj)=doocsread(addr_charge);

    % compared it with last measurement
    if jj > 1
        while ts(jj) == ts(jj-1)
            ddd_read            = doocsread([app.addr_cam, 'IMAGE_EXT_ZMQ']);
            img_sig(:,:,jj)       = ddd_read.data.val_val;
            display( [' - (', num2str(jj), ') same data ... wait ...']);
            pause(1/app.rep_rate)
        end
    end

    pause(1/app.rep_rate)

end


img_mean = squeeze(mean(img_sig, 3));

%%%

% ccharge_7FL2XTDS=,

a=1
length_y, length_x=size(img_sig)
%y_axis, y_axis= collibrataxis(length_y, length_x)
imag=hlc_clean_image( squeeze(img_sig(a:end,:)))

% time
tmp_profile         = img_mean;
[x_com(i), x_var(i), x_fwhm(i), x_axis(i), x_profile(i)]  = get_profile_stats(tmp_profile, timecal_fspixel);

%%
volt=[-1,0, 1]

figure(122321)
plot(volt,x_fwhm, '-.o')

