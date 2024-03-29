% in this fuction we load a polarix immage and the last energy and time
% collibratioan and anaglyse the imgae and save in the output_filename

% Najmeh Mirian, DESY, 29, Nov, 2023

function output_filename= analysing_image(filename, a)
close all
imag_file=load(filename);
load('/home/ttflinac/user/mflorian/PolariX/FL2/time_calibration/time_calib_9FL2XTDS.mat', 'timecal_fspixel', ...
    'timecal_fspixel_err', 'timecal_streak', 'time_res');
load('/home/ttflinac/user/mflorian/PolariX/FL2/energy_calibration/erg_calib_9FL2XTDS.mat', 'ergcal', 'ergcal_err')
display('data is loaded');

addpath('/System/Volumes/Data/home/ttflinac/user/mirian/FL2_funs')

fontSize        = 14;
calib_x             = timecal_fspixel;
calib_y             = ergcal;
timeres_calib       = time_res;


if a <200

img_sig    = imag_file.img(1:end, :, :);

else 
    img_sig    = imag_file.img(a:end, :, :);
end 

num_sig     = size(img_sig, 3);
length_y    = size(img_sig, 1);
length_x    = size(img_sig, 2);



x_com   = zeros([1, num_sig]);   y_com   = zeros([1, num_sig]);
x_var   = zeros([1, num_sig]);   y_var   = zeros([1, num_sig]);
x_fwhm  = zeros([1, num_sig]);   y_fwhm  = zeros([1, num_sig]);
x_profile= zeros([num_sig, length_x]);
y_profile= zeros([num_sig, length_y]);
y=1:length_y;

for i=1:num_sig, charge_7FL2XTDS(i)=imag_file.charge(i).data; end
charge_7FL2XTDS_mean=mean(charge_7FL2XTDS)  ; %V check it
tmp2_img=squeeze(transpose(img_sig(:,:,end)));
%
%%
parfor jj = 1:num_sig

    tmp_img             = hlc_clean_image( squeeze(img_sig(:,:,jj)) );

    img_filt(:,:,jj)    = tmp_img;

    % time
    tmp_profile         = mean(img_filt(:,:,jj));

    [x_com(jj), x_var(jj), x_fwhm(jj), x_axis, x_profile(jj,:)]  = get_profile_stats(tmp_profile, calib_x);

    % erg
    tmp_profile         = mean(img_filt(:,:,jj), 2);
    %     tmp_profile         = medfilt1(tmp_profile, 1);
    [y_com(jj), y_var(jj), y_fwhm(jj), y_axis, y_profile(jj,:)]  = get_profile_stats(tmp_profile', calib_y);

end
tmp2_img=squeeze(transpose(img_filt(:,:,end)));
%% plot 1

erg_axis            = ergcal * ((1:length_y)-length_y/2);
time_axis           = timecal_fspixel * ((1:length_x)-length_x/2);
erg_pos_good =  find(y_profile(end,:)>0); erg_pos_good = min(erg_pos_good):max(erg_pos_good);
time_pos_good = find(x_profile(end,:)>0); time_pos_good = min(time_pos_good):max(time_pos_good);
cent_energy_cr=zeros(num_sig, length_x);
slice_enrgy_spread=zeros(num_sig, length_x);
%%
% energy spread calcaultion
for jj=1:num_sig
    parfor n=min(time_pos_good):max(time_pos_good)


        Ans=hlc_fit_gaussian(erg_pos_good, smoothdata(img_filt(erg_pos_good,n,jj)) );
        %Ans=[baseline, height,mu, sigma]
        if Ans(4)>100
        else

            slice_enrgy_spread(jj, n)=Ans(4)*calib_y;
            cent_energy_cr(jj, n)=Ans(3)*calib_y;
        end
%         if slice_energy_spread(jj, n)> 100
%                         slice_energy_spread(jj, n)=0
%                     end

    end
end


%% jetting of stabilitty correction   ********************NOT applicable wait
for jj=1: num_sig
    ftemp=find(x_profile(jj,:)==max(x_profile(jj,:)));
    mass_center(jj)=ftemp(1);

end



%%
figure(10012)
set(gcf, 'Position', [800, 1, 1500, 1600])

subplot(2,2,1)
imagesc(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'),...
    erg_axis(erg_pos_good), flipud(tmp2_img(time_pos_good,erg_pos_good(end:-1:1))') )
grid on
title(['single image '], 'FontSize', fontSize, 'Interpreter', 'none')
xlabel('time t /fs)', 'FontSize', fontSize)
ylabel('rel. energy deviation \delta / 10^{-3}' , 'FontSize', fontSize)
set(gca, 'FontSize', fontSize,'YDir', 'normal')
lim_y = get(gca, 'YLim'); lim_x = get(gca, 'XLim');

subplot(2,2,2)
p1 = plot( time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'),...
    cent_energy_cr(end,time_pos_good)-mean(cent_energy_cr(end,time_pos_good)));
grid on
title('energy deviation' ,'FontSize', fontSize)
xlabel('norm. erg density', 'FontSize', fontSize)
ylabel('central of energy \delta / 10^{-3}', 'FontSize', fontSize)
legend([p1(1)], {['rms: ' num2str(y_var(end), '%5.2f'), ' \cdot 10^{-3}', 10, ...
    'FWHM: ', num2str(y_fwhm(end), '%5.2f'), ' \cdot 10^{-3}']}, ...
    'Location','northoutside','FontSize', fontSize-2)
set(gca, 'FontSize', fontSize, 'YDir', 'normal');%, 'Ylim', -lim_y)

subplot(2,2,3)

yyaxis left
p1 = plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'),...
    1e-3*charge_7FL2XTDS_mean*1e-9*1e15 * x_profile(end,time_pos_good));
grid on
title('current and energy spread', 'Interpreter', 'none')
xlabel('time t /fs', 'FontSize', fontSize)
ylabel('peak current /kA', 'FontSize', fontSize)
%

legend([p1], {['rms: ' num2str(x_var(end), '%5.2f'), ' fs', 10, ...
    'FWHM: ', num2str(x_fwhm(end), '%5.2f'), ' fs', 10, ...
    'Time res: ' num2str((timeres_calib), '%5.1f'), ' fs']}, ...
    'Location','northoutside','FontSize', fontSize-2)
yyaxis right

p2 = plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'), slice_enrgy_spread(end,time_pos_good))


%slice_enrgy_spread(:,time_pos_good)
ylabel('energy spread (Mev)')
legend([p2(1),p1(1)], { ['energy spread'],['current   ','rms: ' num2str(x_var(end), '%5.2f'), ' fs', 10, ...
    'FWHM: ', num2str(x_fwhm(end), '%5.2f'), ' fs', 10, ...
    'Time res: ' num2str((timeres_calib), '%5.1f'), ' fs']}, ...
    'FontSize', fontSize-2, 'Location','northoutside')


set(gca, 'FontSize', fontSize, 'Xlim', lim_x)
%
%

subplot(2,2,4)
MarkerSizz = 12;
color1 = [0,0.4470, 0.7410];
color2 = [0.85, 0.325, 0.098];
yyaxis left
p1 = plot(1:num_sig, x_var, 'kx', 'MarkerSize', MarkerSizz);
hold on
p2 = plot(1:num_sig, x_fwhm, 'x', 'Color', color1, 'MarkerSize', MarkerSizz);
hold off
ylabel('bunch length \sigma_t /fs', 'FontSize', fontSize)
yyaxis right
p3 = plot(1:num_sig, y_var, 'x', 'Color', color2, 'MarkerSize', MarkerSizz);
ylabel('rms energy spread \sigma_E /10^{-3}', 'FontSize', fontSize)
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = color2;
grid on
xlabel('#measurement', 'FontSize', fontSize)

legend([p1(1), p2(1), p3(1)], {['rms bu length, mean: ', num2str(mean(x_var), '%5.2f'), '\pm', num2str(std(x_var), '%5.2f') ' fs'], ...
    ['fwhm bu length, mean: ', num2str(mean(x_fwhm), '%5.2f'), '\pm', num2str(std(x_fwhm), '%5.2f') ' fs'], ...
    ['rms E spread, mean: ', num2str(mean(y_var), '%5.2f'), '\pm', num2str(std(y_var), '%5.2f') ' \cdot 10^{-3}']}, ...
    'Location','northoutside', 'FontSize', fontSize-2)

set(gca, 'FontSize', fontSize)
%% MBI 
% if num_sig==1
%    subplot(2,2,4)
% FFT2d=abs(fftshift( fft2(flipud(tmp2_img(time_pos_good,erg_pos_good)'))));
% h=pcolor(FFT2d);
% set(h, 'EdgeColor', 'none');
% end 
t=date;
Year=year(t);
fulder_add=['/home/ttflinac/data/PolariX/', num2str(Year), '/'];
timestamp       = datestr(clock, 'yyyy-mm-ddTHHMMSS');
current=1e-3*charge_7FL2XTDS_mean*1e-9*1e15 * x_profile;
save([fulder_add, timestamp, '.mat'],'timestamp' ,  'imag_file',  'slice_enrgy_spread',...
    'cent_energy_cr', 'x_var', 'x_fwhm', 'y_var','ergcal', 'timecal_fspixel','time_pos_good',...
    'current')
output_filename=[fulder_add, timestamp, '.mat']
%% chirp 
%% chirp 
%im = imshow(image);
try 
 clear  time_de cent_energy
 roi = drawrectangle(subplot(2,2,1));
l(1)=roi.Position(1)
 l(2)=roi.Position(3)
             
%index1=find(time_axis==l(1))
[val,index1]=min(abs(time_axis-l(1)));
[val,index2]=min(abs(time_axis-l(2)));
i=1;
for n=index1 : index2

   cent_energy(:,i)= mean(cent_energy_cr(:, n));
   time_de(i)=time_axis(i);
   
   i=i+1;
   
end 

p = polyfit( time_de,cent_energy,3);
y1 = polyval(p,time_de);


figure ()
p1=plot(time_de, cent_energy,'o', 'Color', 'm', 'MarkerSize', 2)
hold on
p2=plot(time_de,y1, 'LineWidth', 2)
legend([p1, p2], { 'data', ['fitting data',10, ...
    ' a=',num2str(p(1)), ' Mev/fs^2', 10  , ...
'b=' ,num2str(p(2)),  'Mev/fs', 10 ...
, ' c=', num2str(p(3)), 'Mev']})
hold off

end 


