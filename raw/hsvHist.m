%% Code for plotting the HSV histogram per person
% Author: Athira Nambiar
% Date 19/12/2013

%%
clear all;
close all;
clc;

%%
pwd;
cd 1/1_1
BCK=imread('rgbframe9999.png')
figure,
subplot(2,3,1),imshow('rgbframe280.png');title('RGB image');       % show the RGB image
subplot(2,3,2),imshow('frame280.png');   title('Depth image');     % show the Depth image
RGB=imread('rgbframe280.png');                                     % read the RGB image
DPT=imread('frame280.png');                                     % read the RGB image
HSV=rgb2hsv(RGB);                                                  % convert RGB to HSV image
subplot(2,3,3),imshow(HSV);title('HSV image');                     % show the HSV image

h = HSV(:, :, 1);                                % Hue image.
s = HSV(:, :, 2);                                % Saturation image.
v = HSV(:, :, 3);                                % Value (intensity) image.

subplot(2,3,4),imshow(h),title('Hue image');
subplot(2,3,5),imshow(s),title('Saturation image');
subplot(2,3,6),imshow(v),title('Value image');

%% Creating the binary mask for the person
%[final Z verticesX verticesY]=mask (RGB,DPT,BCK,width)
