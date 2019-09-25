
function [outim]=chromaticity_(testrgbFrame,backgroundFrame,verticesXX,verticesYY)

% close all;
% clear all;
% clc

% Read the image frames   

% cd 'C:\Users\athira\Desktop\Kinect_Softwares\knect_matlab\21';
% 
% testrgbFrame='rgbframe270.png';
% backgroundFrame='rgbframe340.png';

a=imread(testrgbFrame);       % frame including human body
b=imread(backgroundFrame);     % background image

% figure,
% set(gcf, 'Position', get(0,'Screensize'));

% subplot(1,2,1),imshow(a)
 R_a = a(:,:,1);
 G_a = a(:,:,2);
 B_a = a(:,:,3);
 S_a = sum(a,3);
 Rn_a=double(R_a)./S_a;
 Gn_a=double(G_a)./S_a;
 Bn_a=double(B_a)./S_a;
 img1 = cat(3,Rn_a,Gn_a,Bn_a);
%  subplot(3,2,2),imshow(img1) 
 
% subplot(3,2,3),imshow(b)
 R_b = b(:,:,1);
 G_b = b(:,:,2);
 B_b = b(:,:,3);
 S_b = sum(b,3);
 Rn_b=double(R_b)./S_b;
 Gn_b=double(G_b)./S_b;
 Bn_b=double(B_b)./S_b;

img2 = cat(3,Rn_b,Gn_b,Bn_b);
% subplot(3,2,4),imshow(img2);

diff=20.*(img2-img1);% subplot(3,2,5),imshow(diff)

BW=im2bw(diff);

se = strel('diamond',2);
bw2 = imdilate(BW,se);           % Binary dilation is carried out on thresholded image  
BWfill = imfill(bw2, 'holes');  % Fill the image regions and holes inside
BWnobord = imclearborder(BWfill, 1);  % Clears the image border
% plot(verticesXX, verticesYY, 'm-', 'LineWidth', 2)
% 
outim = bwlargestblob(BWnobord,8);
outim=outim>0;
% subplot(3,2,6),
% imshow(outim);
  
end
