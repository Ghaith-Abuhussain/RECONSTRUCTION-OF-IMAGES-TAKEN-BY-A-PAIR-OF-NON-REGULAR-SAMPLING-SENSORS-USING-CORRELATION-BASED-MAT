clear all ; clc ; close all ;
SENSOR_Left=imread('multi_view/books/left_Books.png');
SENSOR_Right=imread('multi_view/books/right_Books.png');
subplot(1,2,1);imshow(rgb2gray(SENSOR_Left));
subplot(1,2,2);imshow(rgb2gray(SENSOR_Right));
%%
figure
S=stereoAnaglyph(SENSOR_Right,SENSOR_Left);
imshow(S);
title('Red-cyan composite view of the stereo images');

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the parameters
    fse_parameters = struct();
    fse_parameters.block_size = 16;
    fse_parameters.conc_weighting = 0.5;
    fse_parameters.debug = 0;
    fse_parameters.rhos = [0.80, 0.70, 0.66, 0.64];
    fse_parameters.block_size_min = 4;
    fse_parameters.fft_size = 32;
    fse_parameters.max_iter = 100;
    fse_parameters.min_iter = 20;
    fse_parameters.iter_const = 1000;
    fse_parameters.orthogonality_correction = 0.5;
%%
% NON-REGULAR SAMPLING FOR THE LEFT AND THE RIGHT SENSORS
B1=SENSOR_Right;
B2=SENSOR_Left;
[sampled_image1,error_mask1]=non_regular_sampling(B1,4);
[sampled_image2,error_mask2]=non_regular_sampling(B2,4);

%%
%PLOT THE NON-REGULAR SAMPLED IMAGES
subplot(1,2,1);imshow(sampled_image1);title('Right Sensor');
subplot(1,2,2);imshow(sampled_image2);title('Left Sensor');
valid_pixels_Before=find(sampled_image1(:,:,1)~=0);
%%
% EXTRAPOLATE THE LEFT AND RIGHT SAMPLED IMAGES
image_YCbCr = rgb2ycbcr(sampled_image1) .*(error_mask1/255);
    tic();
        reconstructed_image1_YCbCr(:,:,1) = processing_whole_image(image_YCbCr(:,:,1), error_mask1/255, fse_parameters);
        reconstructed_image1_YCbCr(:,:,2) = processing_whole_image(image_YCbCr(:,:,2), error_mask1/255, fse_parameters);
        reconstructed_image1_YCbCr(:,:,3) = processing_whole_image(image_YCbCr(:,:,3), error_mask1/255, fse_parameters);
        reconstructed_img1 = ycbcr2rgb(reconstructed_image1_YCbCr/255);
    toc();
    imwrite(reconstructed_img1,'r1.png');
%
image_YCbCr = rgb2ycbcr(sampled_image2) .*(error_mask2/255);
    tic();
        reconstructed_image2_YCbCr(:,:,1) = processing_whole_image(image_YCbCr(:,:,1), error_mask2/255, fse_parameters);
        reconstructed_image2_YCbCr(:,:,2) = processing_whole_image(image_YCbCr(:,:,2), error_mask2/255, fse_parameters);
        reconstructed_image2_YCbCr(:,:,3) = processing_whole_image(image_YCbCr(:,:,3), error_mask2/255, fse_parameters);
        reconstructed_img2 = ycbcr2rgb(reconstructed_image2_YCbCr/255);
    toc();
     imwrite(reconstructed_img2,'r2.png');
%%
% PLOT THE RESULTING IMAGES
subplot(1,2,1);imshow(reconstructed_img1);title('Right Before');
subplot(1,2,2);imshow(reconstructed_img2);title('Left Before');

%%
% FIND THE DISPARITY MAP LEFT TO RIGHT AND RIGHT TO LEFT
tic();
y=disparity_map(reconstructed_img2,reconstructed_img1,64,1,0.9);
y1=disparity_map(reconstructed_img1,reconstructed_img2,64,2,0.9);
toc();

D_RightToLeft=y1;
D_LeftToRight=y;

subplot(2,2,1);imshow(rgb2gray(SENSOR_Left));title('Left');
subplot(2,2,2);imshow(rgb2gray(SENSOR_Right));title('Right');
subplot(2,2,3);imshow(mat2gray(D_RightToLeft));title('Right To Left');
subplot(2,2,4);imshow(mat2gray(D_LeftToRight));title('Left To Right');
%%
% PROJECT THE VALID PIXELS LEFT TO RIGHT AND RIGHT TO LEFT DEPENDING ON
% THE DISPARITY MAPS
sampled_image_Left=sampled_image2;
sampled_image_Right=sampled_image1;
error_maskRight=error_mask1;
error_maskLeft=error_mask2;
[M1,N1]=size(sampled_image2(:,:,1));
for i=1:M1
    for j=1:N1
        if sampled_image2(i,j,:) ~=0
            if j-D_LeftToRight(i,j) >=1 && D_LeftToRight(i,j)>0 && j-D_LeftToRight(i,j)<=N1
                if sampled_image2(i,j-D_LeftToRight(i,j),:)==0;
                    sampled_image_Right(i,j-D_LeftToRight(i,j),:) = sampled_image_Left(i,j,:);
                    error_maskRight(i,j-D_LeftToRight(i,j),:)=255;
                end
            end
        end
        if sampled_image1(i,j,:) ~=0
            if j+D_LeftToRight(i,j) <= N1 && D_LeftToRight(i,j)>0 && j+D_LeftToRight(i,j) >= 1
                if sampled_image1(i,j+D_LeftToRight(i,j),:)==0;
                    sampled_image_Left(i,j+D_LeftToRight(i,j),:) = sampled_image_Right(i,j,:);
                    error_maskLeft(i,j+D_LeftToRight(i,j),:)=255;
                end
            end
        end
    end
end
%%
% SHOW THE NEW SAMPLED IMAGES
subplot(1,2,1);imshow(sampled_image_Right);title('Right Sensor');
subplot(1,2,2);imshow(sampled_image_Left);title('Left Sensor');
valid_pixels_After=find(sampled_image_Right(:,:,1)~=0);
%%
% EXTRAPOLATE THE NEW SAMPLED IMAGES
image_YCbCr = rgb2ycbcr(sampled_image_Right) .*(error_maskRight/255);
    tic();
        reconstructed_imageRight_YCbCr(:,:,1) = processing_whole_image(image_YCbCr(:,:,1), error_maskRight/255, fse_parameters);
        reconstructed_imageRight_YCbCr(:,:,2) = processing_whole_image(image_YCbCr(:,:,2), error_maskRight/255, fse_parameters);
        reconstructed_imageRight_YCbCr(:,:,3) = processing_whole_image(image_YCbCr(:,:,3), error_maskRight/255, fse_parameters);
        reconstructed_imgRight = ycbcr2rgb(reconstructed_imageRight_YCbCr/255);
    toc();
    
image_YCbCr = rgb2ycbcr(sampled_image_Left) .*(error_maskLeft/255);
    tic();
        reconstructed_imageLeft_YCbCr(:,:,1) = processing_whole_image(image_YCbCr(:,:,1), error_maskLeft/255, fse_parameters);
        reconstructed_imageLeft_YCbCr(:,:,2) = processing_whole_image(image_YCbCr(:,:,2), error_maskLeft/255, fse_parameters);
        reconstructed_imageLeft_YCbCr(:,:,3) = processing_whole_image(image_YCbCr(:,:,3), error_maskLeft/255, fse_parameters);
        reconstructed_imgLeft = ycbcr2rgb(reconstructed_imageLeft_YCbCr/255);
    toc();
%%
% BEFORE AND AFTER FOR RIGHT SENSOR
subplot(1,2,1);imshow(reconstructed_imgRight);title('Right After');
subplot(1,2,2);imshow(reconstructed_img1);title('Right Before');
%%
% BEFORE AND AFTER FOR LEFT SENSOR
subplot(1,2,1);imshow(reconstructed_imgLeft);title('Left After');
subplot(1,2,2);imshow(reconstructed_img2);title('Left Before');
%%
fprintf('PSNR:  %0.4f dB\n', psnr(im2double(C), reconstructed_img));
fprintf('SSIM:   %0.4f\n', ssim(im2double(C), reconstructed_img));







