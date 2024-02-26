function get_off_lasing_data(filename,a )

imag_lasing_off=load(filename);
% calibration data
load('/home/ttflinac/user/mflorian/PolariX/FL2/time_calibration/time_calib_9FL2XTDS.mat', 'timecal_fspixel', ...
    'timecal_fspixel_err', 'timecal_streak', 'time_res');
load('/home/ttflinac/user/mflorian/PolariX/FL2/energy_calibration/erg_calib_9FL2XTDS.mat', 'ergcal', 'ergcal_err')
display('data is loaded');

%num_sig=app.num_shot

%img_sig    = imag_file.img(end:a, :, :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% analysis
%img_bgr_mean=Bakg.img_bgr_mean;
% img_sig=imag_lasing_on.img;
% charge_7FL2XTDS_mean=mean(mag_lasing_on.charge.data)  ; %V check it
fontSize        = 14;
calib_x             = timecal_fspixel;


calib_y             = ergcal;

timeres_calib       = time_res;
%%
img_sig=imag_lasing_off.img(end:a, :, :);
num_sig     = size(img_sig, 3);
for i=1:num_sig, charge_7FL2XTDS(i)=imag_lasing_off.charge(i).data; end

charge_7FL2XTDS_mean=mean(charge_7FL2XTDS)  ; %V check it


img_filt    = img_sig;

num_sig     = size(img_sig, 3);
length_y    = size(img_sig, 1);
length_x    = size(img_sig, 2);

%img_bgr_mean = squeeze(mean(img_bgr, 3));

x_com       = zeros([1, num_sig]);              y_com       = zeros([1, num_sig]);
x_var       = zeros([1, num_sig]);              y_var       = zeros([1, num_sig]);
x_fwhm      = zeros([1, num_sig]);              y_fwhm      = zeros([1, num_sig]);
x_profile   = zeros([num_sig, length_x]);       y_profile   = zeros([num_sig, length_y]);
y=1:length_y;
%%

parfor jj = 1:num_sig

    tmp_img             = hlc_clean_image( squeeze(img_sig(:,:,jj)) );
    %tmp_img             = get_ROI(squeeze(img_sig(:,:,jj)) - img_bgr_mean);
    img_filt(:,:,jj)    = medfilt2(tmp_img);

    % time
    tmp_profile         = mean(img_filt(:,:,jj));
    %     tmp_profile(2010)=tmp_profile(2009); % remove hot pixel
    %     tmp_profile         = medfilt1(tmp_profile, 1);
    [x_com(jj), x_var(jj), x_fwhm(jj), x_axis, x_profile(jj,:)]  = get_profile_stats(tmp_profile, calib_x);

    % erg
    tmp_profile         = mean(tmp_img, 2);
    %     tmp_profile         = medfilt1(tmp_profile, 1);
    [y_com(jj), y_var(jj), y_fwhm(jj), y_axis, y_profile(jj,:)]  = get_profile_stats(tmp_profile', calib_y);

end
current=1e-3*charge_7FL2XTDS_mean*1e-9*1e15*x_profile(jj,time_pos_good)

tmp2_img=squeeze(transpose(img_filt(:,:,end)));

%% plot 1

erg_axis            = ergcal * ((1:length_y)-length_y/2);
time_axis           = timecal_fspixel * ((1:length_x)-length_x/2);
erg_pos_good =  find(y_profile(end,:)>0); erg_pos_good = min(erg_pos_good):max(erg_pos_good);
time_pos_good = find(x_profile(end,:)>0); time_pos_good = min(time_pos_good):max(time_pos_good);
cent_energy_cr=zeros(num_sig, length_x);
slice_enrgy_spread=zeros(num_sig, length_x);
%%
%     figure; plot(slice_enrgy_spread_ave)
%%
for jj=1:num_sig
    for n=min(time_pos_good):max(time_pos_good)

        %  mu0=  find( img_filt(erg_pos_good,n,jj)==max(img_filt(erg_pos_good,n,jj)));
        %  [sigma, mu] = gaussfit( erg_pos_good, smoothdata(img_filt(erg_pos_good,n,jj)) , 0, erg_pos_good(mu0(1)));
        Ans=hlc_fit_gaussian(erg_pos_good, smoothdata(img_filt(erg_pos_good,n,jj)) );
        %Ans=[baseline, height,mu, sigma]
        if isnan(Ans(3))
        else

            slice_enrgy_spread(jj, n)=Ans(4)*calib_y;
            cent_energy_cr(jj, n)=Ans(3)*calib_y;
        end

    end
end
