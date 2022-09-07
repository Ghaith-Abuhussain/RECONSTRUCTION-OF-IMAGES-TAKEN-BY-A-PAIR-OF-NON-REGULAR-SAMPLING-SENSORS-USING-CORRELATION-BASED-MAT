% READ THE IMAGES
clear all ; clc ; close all ;
SENSOR_Left=imread('multi_view/books/left_Books.png');
SENSOR_Right=imread('multi_view/books/right_Books.png');
subplot(1,2,1);imshow(rgb2gray(SENSOR_Left));
subplot(1,2,2);imshow(rgb2gray(SENSOR_Right));
%%
% FIND DISPARITY MAPS
tic();
y=disparity_map(SENSOR_Left,SENSOR_Right,4,1);
y1=disparity_map(SENSOR_Right,SENSOR_Left,4,2);
toc();
%%
%PLOT THE RESULT
subplot(2,2,1);imshow(SENSOR_Left);title('SENSOR LEFT');
subplot(2,2,2);imshow(SENSOR_Right);title('SENSOR RIGHT');
subplot(2,2,3);imshow(mat2gray(y));title('disparity map RightToLeft');
subplot(2,2,4);imshow(mat2gray(y1));title('disparity map LeftToRight');