function clean_img = hlc_clean_image(input_img)
% Clean the background from the given image.
%
% clean_img = hlc_clean_image(input_img)
%   This function cleans the given image from isolated bright pixels and
%   automatically masks out all portions of the image that the algorithm
%   identifies as background noise.
%   The algorithm is based on an algorithm by Luca Giannessi.
%
% -- Lars Froehlich, Luca Giannessi

    input_img = hlc_denoise_image(input_img);

    gaussian_filter = calc_gaussian_2d(35, 7);
    img_filtered = filter2(gaussian_filter, input_img);

    good_soglia = ricerca_soglia_lg(img_filtered, 0.1, 0.03, 35);

    [~, maschera] = mask_image(img_filtered, good_soglia);
    fondo_punto = min(min(input_img));
    input_img = input_img - fondo_punto;

    clean_img = double(input_img) .* maschera;
end

function matrix = calc_gaussian_2d(pixels, sigma)
    gauss1d = hlc_calc_gaussian(0:pixels-1, [0, 1, (pixels-1)/2, sigma]);
    matrix = gauss1d' * gauss1d; 
end

function threshold = ricerca_soglia_lg(img_input, guess, tol, n_iter)
    area_s = size(img_input, 1) * size(img_input, 2);
    signal_s = sum(sum(img_input));

    for i = 1:n_iter
        [img_k, mask_k] = mask_image(img_input, guess);

        area_k = sum(sum(mask_k));
        signal_k = sum(sum(img_k));

        noise_k = (signal_s - signal_k) / (area_s - area_k + 1);
        signal_d = signal_k/(area_k + 1);
        s_n_ratio = abs(signal_d - noise_k) / (noise_k + 1);
        new_guess = (1 + 1/s_n_ratio) / s_n_ratio / sqrt(2);
        
        if abs(new_guess - guess) < tol * guess
            threshold = min([new_guess, 0.3]);
            return;
        else
            guess = new_guess;
        end
    end
    
    threshold = guess;
end

function [masked_image, mask] = mask_image(img_in, fraction_of_maximum)
    max_intensity = max(max(img_in));
    mask = (img_in >= max_intensity * fraction_of_maximum);
    masked_image = img_in .* mask;
end
