clear all ; clc ; close all ; 
A=imread('img/400.jpg');
% the parameters
    fsr_parameters = struct();
    
    % default variables are set in this function.
    fsr_parameters.block_size = 16;
    fsr_parameters.conc_weighting = 0.5;
    fsr_parameters.debug = 0;
    fsr_parameters.rhos = [0.80, 0.70, 0.66, 0.64];
    fsr_parameters.block_size_min = 4;
    fsr_parameters.fft_size = 32;
    fsr_parameters.max_iter = 100;
    fsr_parameters.min_iter = 20;
    fsr_parameters.iter_const = 1000;
    fsr_parameters.orthogonality_correction = 0.5;
    [sampled_image1,error_mask1]=non_regular_sampling(A,4);
%%
    subplot(1,2,1);imshow(A);title('the LR image');
    subplot(1,2,2);imshow(sampled_image1);title('non regular sampled image');
%%
    image_YCbCr = rgb2ycbcr(sampled_image1) .*(error_mask1/255);
    tic();
        reconstructed_image1_YCbCr(:,:,1) = processing_whole_image(image_YCbCr(:,:,1), error_mask1/255, fsr_parameters);
        reconstructed_image1_YCbCr(:,:,2) = processing_whole_image(image_YCbCr(:,:,2), error_mask1/255, fsr_parameters);
        reconstructed_image1_YCbCr(:,:,3) = processing_whole_image(image_YCbCr(:,:,3), error_mask1/255, fsr_parameters);
        reconstructed_img1 = ycbcr2rgb(reconstructed_image1_YCbCr/255);
    toc();
    %%
    % results :
    C = imread('img/800.jpg');
    subplot(1,3,1);imshow(C);title('reference image');
    subplot(1,3,3);imshow(reconstructed_img1);title('reconstructed image');
    subplot(1,3,2);imshow(sampled_image1);title('non regular sampled image');

    fprintf('PSNR:  %0.4f dB\n', psnr(im2double(C),reconstructed_img1));
    fprintf('SSIM:   %0.4f\n', ssim(im2double(C),reconstructed_img1));