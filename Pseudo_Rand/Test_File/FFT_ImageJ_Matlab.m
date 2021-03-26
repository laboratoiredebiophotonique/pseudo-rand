clear;clc;

Im = imread('circ.png');
Im = rgb2gray(Im);
figure, imshow(Im)
A = fft2(double(Im)/255);
Ashifted = fftshift(A);

Ashifted(97,158)
Ashifted(93,165)