%%Code for background subtraction and the extration of only the foreground
% of the new image- in our case, get the RGB image of the person.
% Author: Athira
%Date: 16-12-2013

clc;        % Clear command window.
clear all;      % Delete all variables.
close all;	% Close all figure windows .

workspace;	% Make sure the workspace panel is showing.
fontSize = 16;  
tole=25;
width=15;

HSV_person=[];
addpath(genpath('C:\HSV_MoG\final codes'));

%% !!!!!!!!!!!!!!! DATA AQUISITION !!!!!!!!!!!!!!

% Read the image frames   

cd 'C:\HSV_MoG\MoG\DB4\1\1'

rgb_imagefiles = dir('rgbframe*.png');                  % get list of rgb .png files in this directory
noOfRgbFrames = length(rgb_imagefiles);                 % Number of RGB files found
depth_imagefiles = dir('frame*.png');                   % get list of depth .png files in this directory
noOfDepthFrames = length(depth_imagefiles);             % Number of depth files found
backgroundFrame='rgbframe9999.png';

for k =40:60
%k=86 (noOfRgbFrames-1)
testdepthFrame=depth_imagefiles(k).name;                %Read the k'th depth frame
testrgbFrame=rgb_imagefiles(k).name;                    %Read the k'th rgb frame

a=imread(testrgbFrame);                                 % frame including human body
b=imread(backgroundFrame);                              % background image

a=rgb2gray(a);                                          % Converting to gray images
b=rgb2gray(b);

a=imadjust(a, stretchlim(a), [0 1]);                    % Illumination adjustment
b=imadjust(b, stretchlim(b), [0 1]);

diff=b-a;                                               % Difference of RGB images

        figure,
        subplot(3,3,1)
        imshow(diff); title('difference of RGB image')% difference of RGB images

% % Grayscale and BW conversion

im=diff   ;                                             % grayscale of difference of RGB image
level = graythresh(im);                                 % calculate the smart threshold level using Otsu's method 
d = im2bw(im,level)   ;                                 % RGB difference image is thresholded and BW image is plotted
          subplot(3,3,2)
          imshow(d), title('Thresholded BW image'),hold on,

% Morphological operations

se = strel('diamond',8);
bw2 = imdilate(d,se);                                   % Binary dilation is carried out on thresholded image  
BWfill = imfill(bw2, 'holes');                          % Fill the image regions and holes inside
BWnobord = imclearborder(BWfill, 1);                    % Clears the image border
          subplot(3,3,3)
          imshow(BWnobord), title('Dilated'),hold on,

%Label each blob so we can make measurements of it
[labeledImage numberOfBlobs] = bwlabel(BWnobord, 8);
% Get all the blob properties.
blobMeasurements = regionprops(labeledImage, 'BoundingBox','Area');
allBlobAreas = [blobMeasurements.Area];
for k = 1 :numberOfBlobs
boundingBox = blobMeasurements(k).BoundingBox;	 % Get box.
aspectRatio(k) = boundingBox(3) / boundingBox(4);
% fprintf('For blob #%d, area = %d, aspect ratio = %.2f\n', ...
% k, allBlobAreas(k), aspectRatio(k));
end

% Find the biggest binary blob and plot the bounding box
[r,c] = find(allBlobAreas==max(allBlobAreas(:))); 
s=blobMeasurements(c).BoundingBox;
x1 = s(1);
y1 = s(2);
x2 = x1 + s(3) - 1;
y2 = y1 + s(4) - 1;
verticesX = [x1 x2 x2 x1 x1];
verticesY = [y1 y1 y2 y2 y1];
% plot(verticesX, verticesY, 'r-', 'LineWidth', 2);


% Generate the depth image

i=imread(testdepthFrame);
I=i(:,:,1);
        subplot(3,3,4)
        imshow(I);hold on;
        plot(verticesX, verticesY, 'r-', 'LineWidth', 2);
        title('DepthImage')

% generate the adapting shifing for the bounding box in order to cover the
% human body in the depth image ...!
% med=mode(mode(double(I(y1:y2,x1:x2))))
human_box=((double(I(y1:y2,x1:x2)))).*BWnobord(y1:y2,x1:x2);
human_point=mode(human_box(human_box~=0));
med=human_point;


max_threshold=med+width;
min_threshold=med-width;
human = I; 
human = human<max_threshold & human>min_threshold;
human(:,1:x1-tole)=0;
human(:,x2+tole:end)=0;
human(1:y1-3*tole,:)=0;
human(y2+tole:end,:)=0;
human_depth=(double(I).*human);
            subplot(3,3,5)
            imshow(human_depth);hold on;
%             plot(verticesX_new, verticesY_new, 'm-', 'LineWidth', 2);
            title('Depth silhouette ')
Z=human_depth;

        subplot(3,3,6)
        imshow(testrgbFrame),title('rgb image')

%Figure out the final human silhouette by combining both the depth
%silhouette with the foreground silhouette
        subplot(3,3,7)
        final=Z.*BWnobord;
        imshow(final)
        title(' human silhouette extracted in RGB colorspace')
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mask in chromaticity space also and combine with RGB mask to obtain the final composite mask

%final=final>0;                              % Final is the output from the function bsub- human silhouette in RGB space
final(:,1:verticesX(1))=0;
final(:,verticesX(2):end)=0;
final(1:verticesY(1),:)=0;
final(verticesY(3):end,:)=0;

% subplot(1,3,1)
% imshow(final),hold on,
% plot(verticesX, verticesY, 'r-', 'LineWidth', 2);
% plot(X0,Y0,'-m*');
% title(' Human mask generated by background subtracton')

mid=round((((verticesY(3)-verticesY(1))*2/3))+verticesY(1));
weight1=ones(size(final));
weight1(mid+1:end,:)=0;

final=final.*weight1;

[outim]=chromaticity(testrgbFrame,backgroundFrame);  % human silhouette in chromaticity space
outim(:,1:verticesX(1))=0;
outim(:,verticesX(2):end)=0;
outim(1:verticesY(1),:)=0;
outim(verticesY(3):end,:)=0;

% imshow(outim),hold on,
% plot(X0,Y0,'-m*');
% plot(verticesX, verticesY, 'r-', 'LineWidth', 2);
% title(' Human chromaticity mask ')
% 
weight2=ones(size(outim));
weight2(1:mid-1,:)=0;

outim=outim.*weight2;   % Chromaticity mask
        subplot(3,3,8)
        imshow(outim),title('chromaticity mask')

silhouette=final+outim;    % combine the masks in RGB colorspace and chromaticity space to get the composite mask
        subplot(3,3,9)
        imshow(silhouette), title('final mask')
        
%         subplot(2,3,4)
%         imshow(testrgbFrame),title('rgb image')

silhouette(find(silhouette~=0))=1;
%         subplot(2,3,5)
%         imshow(silhouette), title('binary mask')

% Project the binary mask over the RGB image

a=imread(testrgbFrame);       % frame including human body

a_new = a.*repmat(uint8(silhouette),[1,1,3]);
%         subplot(2,3,6), 
%         imshow(a_new),title('foreground extracted')

%reference pages:
%http://www.mathworks.com/matlabcentral/answers/38547
%http://www.mathworks.com/matlabcentral/fileexchange/28512-simple-color-detection-by-hue/content/SimpleColorDetectionByHue.m


%% Convert to HSV image 
HSV=rgb2hsv(a_new);  % convert RGB to HSV image
        figure;
        subplot(2,3,1),imshow(HSV);title('HSV image');                     % show the HSV image

h = HSV(:, :, 1);                                % Hue image.
s = HSV(:, :, 2);                                % Saturation image.
v = HSV(:, :, 3);                                % Value (intensity) image.

        subplot(2,3,2),imshow(h),title('Hue image');
        subplot(2,3,3),imshow(s),title('Saturation image');
        subplot(2,3,4),imshow(v),title('Value image');

numberOfBins_h = 16; 
numberOfBins_s = 16;
numberOfBins_v = 4;
% Get the histogram of the H channel.
[countsH valuesH] = hist(find(h(:)~=0), numberOfBins_h);
%[countsH valuesH] = imhist(h, numberOfBins_h);
% countsH(1)=0
        figure;
        subplot(2, 2,1);
        bar(valuesH, countsH, 'r');
        grid on;
        xlabel('Hue Value'); 
        ylabel('Pixel Count'); 
        title('Histogram of the H Channel', 'FontSize', 15);	
[countsS, valuesS] = hist(find(s(:)~=0), numberOfBins_s);
%[countsS, valuesS] = hist(s, numberOfBins_s);
% countsS(1)=0
        subplot(2, 2, 2);
       bar(valuesS, countsS, 'g');
        xlabel('Saturation Value'); 
        ylabel('Pixel Count'); 
        title('Histogram of the S Channel', 'FontSize', 15);
% Get the histogram of the V channel.
 [countsV, valuesV] = hist(find(v(:)~=0), numberOfBins_v);
%[countsV, valuesV] = hist(v, numberOfBins_v);
% countsV(1)=0
        subplot(2, 2, 3);
        bar(valuesV, countsV, 'b');
        xlabel('Value Value'); 
        ylabel('Pixel Count'); 
        title('Histogram of the V Channel', 'FontSize', 15);

HSV_vector=[countsH,countsS,countsV]
HSV_vector=(HSV_vector/sum(HSV_vector))*1000
HSV_person=[HSV_person;HSV_vector];
 end

save('HSV_person.mat','HSV_person');
%change line 30 for loop and here %'s, line 31 delete, also the figure&
%subplot options comment


