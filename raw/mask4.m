%%Code for background subtraction and the extration of only the foreground
% of the new image- in our case, get the RGB image of the person.
% Author: Athira
%Date: 16-12-2013
% Final version.

clc;        % Clear command window.
clear all;      % Delete all variables.
close all;	% Close all figure windows .

workspace;	% Make sure the workspace panel is showing.
fontSize = 16;  
tole=25;
width=15;

HSV_person=[];

addpath(genpath('C:\Users\AthiraNambiar\Desktop\PHD\codes\2014'));
%% !!!!!!!!!!!!!!! DATA AQUISITION !!!!!!!!!!!!!!

% Read the image frames   

cd  'C:\Users\AthiraNambiar\Desktop\PHD\codes\2014\HSVhistogram\HSVDatabase\8'

rgb_imagefiles = dir('rgbframe*.png');                  % get list of rgb .png files in this directory
noOfRgbFrames = length(rgb_imagefiles)                 % Number of RGB files found
depth_imagefiles = dir('frame*.png');                   % get list of depth .png files in this directory
noOfDepthFrames = length(depth_imagefiles)             % Number of depth files found
backgroundFrame='rgbframe9999.png';

for fr =100%:(noOfRgbFrames-1)
    
testdepthFrame=depth_imagefiles(fr).name;                %Read the k'th depth frame
testrgbFrame=rgb_imagefiles(fr).name;                    %Read the k'th rgb frame

a=imread(testrgbFrame);                                 % frame including human body
b=imread(backgroundFrame);                              % background image

a=rgb2gray(a);                                          % Converting to gray images
b=rgb2gray(b);

a=imadjust(a, stretchlim(a), [0 1]);                    % Illumination adjustment
b=imadjust(b, stretchlim(b), [0 1]);

diff=b-a;                                               % Difference of RGB images

        figure,
        subplot(3,3,1)
        imshow(testrgbFrame),title('RGB image','FontSize', 11)
        
        subplot(3,3,2)
        imshow(diff); title('Difference of RGB image','FontSize', 11)% difference of RGB images

% % Grayscale and BW conversion

im=diff   ;                                             % grayscale of difference of RGB image
level = graythresh(im);                                 % calculate the smart threshold level using Otsu's method 
d = im2bw(im,level)   ;                                 % RGB difference image is thresholded and BW image is plotted
          subplot(3,3,3)
          imshow(d), title('Thresholded BW image','FontSize', 11),hold on,

% Morphological operations

se = strel('diamond',20);
bw2 = imdilate(d,se);                                   % Binary dilation is carried out on thresholded image  
BWfill = imfill(bw2, 'holes');                          % Fill the image regions and holes inside
BWnobord = imclearborder(BWfill, 1);                    % Clears the image border
          subplot(3,3,4)
          imshow(BWnobord), title('Dilated BW image','FontSize', 11),hold on,

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
plot(verticesX, verticesY, 'r-', 'LineWidth', 2);


% Generate the depth image

i=imread(testdepthFrame);
I=i(:,:,1);
        subplot(3,3,5)
        imshow(I);hold on;
        plot(verticesX, verticesY, 'r-', 'LineWidth', 2);
        title('DepthImage','FontSize', 11)

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
            subplot(3,3,6)
            imshow(human_depth);hold on;
            title('Depth silhouette','FontSize', 11)
            Z=human_depth;
    
        

%Figure out the final human silhouette by combining both the depth
%silhouette with the foreground silhouette
        subplot(3,3,7)
        final=Z.*BWnobord;
        imshow(final),hold on
        title(' Human silhouette extracted in RGB colorspace','FontSize', 11)
%         saveas(h,['silhouette' num2str(fr) '.fig'])                                                   
%         close

         %Label each blob so we can make measurements of it
        [labeledImage numberOfBlobs] = bwlabel(final, 8);
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
        xx1 = s(1);
        yy1 = s(2);
        xx2 = xx1 + s(3) - 1;
        yy2 = yy1 + s(4) - 1;
        verticesXX = [xx1 xx2 xx2 xx1 xx1];
        verticesYY = [yy1 yy1 yy2 yy2 yy1];
%         plot(verticesXX, verticesYY, 'm-', 'LineWidth', 2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mask in chromaticity space also and combine with RGB mask to obtain the final composite mask

%final=final>0;                              % Final is the output from the function bsub- human silhouette in RGB space
final(:,1:verticesXX(1))=0;
final(:,verticesXX(2):end)=0;
final(1:verticesYY(1),:)=0;
final(verticesYY(3):end,:)=0;

% subplot(1,3,1)
% imshow(final),hold on,
% plot(verticesX, verticesY, 'r-', 'LineWidth', 2);
% plot(X0,Y0,'-m*');
% title(' Human mask generated by background subtracton')

mid=round((((verticesYY(3)-verticesYY(1))*2/3))+verticesYY(1))
weight1=ones(size(final));
weight1(mid+1:end,:)=0;
                
final=final.*weight1;  % RGB with weight 1 in the upper bogy region

[outim]=chromaticity(testrgbFrame,backgroundFrame,verticesXX,verticesYY);  % human silhouette in chromaticity space
outim(:,1:round(verticesXX(1)))=0;
outim(:,round(verticesXX(2)):end)=0;
outim(1:round(verticesYY(1)),:)=0;
outim(round(verticesYY(3)):end,:)=0;

% imshow(outim),hold on,
% plot(X0,Y0,'-m*');
% plot(verticesX, verticesY, 'r-', 'LineWidth', 2);
% title(' Human chromaticity mask ')

weight2=ones(size(outim));
weight2(1:mid-1,:)=0;

outim2=outim.*weight2;   % Chromaticity mask
        subplot(3,3,8)
        imshow(outim),title('Human silhouette extracted in chromaticity space','FontSize', 11)

silhouette=final+outim2;    % combine the masks in RGB colorspace and chromaticity space to get the composite mask

 %Label each blob so we can make measurements of it
        [labeledImage numberOfBlobs] = bwlabel(silhouette, 8);
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
        xxx1 = s(1);
        yyy1 = s(2);
        xxx2 = xxx1 + s(3) - 1;
        yyy2 = yyy1 + s(4) - 1;
        verticesXXX = [xxx1 xxx2 xxx2 xxx1 xxx1];
        verticesYYY = [yyy1 yyy1 yyy2 yyy2 yyy1];
%         plot(verticesXX, verticesYY, 'g-', 'LineWidth', 2);
        silhouette(:,1:verticesXX(1))=0;
        silhouette(:,verticesXX(2):end)=0;
        silhouette(1:verticesYY(1),:)=0;
        silhouette(verticesYY(3):end,:)=0;

          subplot(3,3,9)
          imshow(silhouette), title('Final mask','FontSize', 11),hold on
%           plot(verticesX, verticesY, 'r-', 'LineWidth', 2);
%           plot(verticesXX, verticesYY, 'm-', 'LineWidth', 2);
%           plot(verticesXXX, verticesYYY, 'g-', 'LineWidth', 2);
          h=figure,imshow(silhouette)
          saveas(h,['a' num2str(fr) '.fig']) 
                   
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
        subplot(1,4,1),imshow(HSV);title('HSV image','FontSize', 15);                     % show the HSV image

h = HSV(:, :, 1);                                % Hue image.
s = HSV(:, :, 2);                                % Saturation image.
v = HSV(:, :, 3);                                % Value (intensity) image.

        subplot(1,4,2),imshow(h),title('Hue image','FontSize', 15);
        subplot(1,4,3),imshow(s),title('Saturation image','FontSize', 15);
        subplot(1,4,4),imshow(v),title('Value image','FontSize', 15);

numberOfBins_h = 16; 
numberOfBins_s = 16;
numberOfBins_v = 4;
% Get the histogram of the H channel.
[countsH valuesH] = hist(h(logical(silhouette)), numberOfBins_h);

% countsH2=countsH/sum(countsH);
%[countsH valuesH] = imhist(h, numberOfBins_h);
% countsH(1)=0
        figure;
        subplot(1,3,1);
        bar(valuesH, countsH, 'r');
        grid on;
        xlabel('Hue Value'); 
        ylabel('Pixel Count'); 
        title('Histogram of the H Channel','FontSize', 15);	
        
 [countsS, valuesS] = hist(s(logical(silhouette)), numberOfBins_s);
% countsS2=countsS/sum(countsS);
%[countsS, valuesS] = hist(s, numberOfBins_s);
% countsS(1)=0
        subplot(1,3, 2);
       bar(valuesS, countsS, 'g');grid on;
        xlabel('Saturation Value'); 
        ylabel('Pixel Count'); 
        title('Histogram of the S Channel','FontSize', 15);
% Get the histogram of the V channel.

 [countsV, valuesV] = hist(v(logical(silhouette)), numberOfBins_v);
%  countsV2=countsV/sum(countsV);
%[countsV, valuesV] = hist(v, numberOfBins_v);
% countsV(1)=0
        subplot(1,3, 3);
        bar(valuesV, countsV, 'b');grid on;
        xlabel('Value Value'); 
        ylabel('Pixel Count'); 
        title('Histogram of the V Channel','FontSize', 15);
        
        
 HSV_vector=[countsH,countsS,countsV];
 HSV_vector=(HSV_vector/sum(HSV_vector)); % Normalize at the end
 HSV_person=[HSV_person;HSV_vector];

%% Trial run to check the differences while normalizing the whole histogram at the end(HSV_vector)
% and each channel separately(HSV_vector2)
% Result: HSV_vector2 is 3 times the HSV_vector

% HSV_vector=[countsH,countsS,countsV];
% HSV_vector2=[countsH2,countsS2,countsV2];
% HSV_vector=(HSV_vector/sum(HSV_vector))*3;
% HSV_vector2=(HSV_vector2)
% figure,
% plot(HSV_vector,'r'); hold on
% plot(HSV_vector2,'b')
%close all;
 end
% 
save('HSV_person.mat','HSV_person');


% %% Code to select the keyframe for shape context Biometrics
% 
% % Save the better silhouette frame (Probably Zero phase) as a .png file.
% 
% I=imcrop
% imshow(I)
% size(I)
% imwrite(I,'SC5_2.png','png')
