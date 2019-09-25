clc,
close all,
clear all, 

scale = 10; M = 40; str = '';x=70;
R = 100; F = 12; N2 = 5; s = 'l';
H1=[];H2=[];H3=[];H4=[];
PTS=[];

X = 255*ones(300,300);
%% Get the avatar silhouette

% PERSON1
I=imread('alexis1_1.png'); %# Load the image data
level = graythresh(I);                 %# Compute an appropriate threshold
BW = im2bw(I,level);                   %# Convert grayscale to binary
[PTS]=contour(BW)                       %# NN approach to get ordered shape contour
my_contour1 = PTS(:,1) +1i*PTS(:,2)
Interp_contour1  = interpolation1(my_contour1,scale,M,str)
histograms1 = extract_features(Interp_contour1,M,R,F,N2,s); % Histograms / Shape Context descriptors of 12x5 bins
    for i=1:M
    hist1=reshape(histograms1(:,:,i),1,[]);; % Reshaped single descriptor of 60D for each shape context
    H1=vertcat(H1,hist1); % 60D distribution for all points
    end
    figure,
    imagesc(~histograms1(:,:,i)');colormap gray


% PERSON2
I=imread('alexis1_3.png'); %# Load the image data
level = graythresh(I);                 %# Compute an appropriate threshold
BW = im2bw(I,level);                   %# Convert grayscale to binary
[PTS]=contour(BW)                       %# NN approach to get ordered shape contour
my_contour2 = PTS(:,1) +1i*PTS(:,2)
Interp_contour2  = interpolation1(my_contour2,scale,M,str)
histograms2 = extract_features(Interp_contour2,M,R,F,N2,s); % Histograms / Shape Context descriptors of 12x5 bins
    for i=1:M
    hist1=reshape(histograms2(:,:,i),1,[]);; % Reshaped single descriptor of 60D for each shape context
    H2=vertcat(H2,hist1); % 60D distribution for all points
    end
    figure,
    imagesc(~histograms2(:,:,i)');colormap gray

% PERSON3
I=imread('alexis8_1.png'); %# Load the image data
level = graythresh(I);                 %# Compute an appropriate threshold
BW = im2bw(I,level);                   %# Convert grayscale to binary
[PTS]=contour(BW)                       %# NN approach to get ordered shape contour
my_contour3 = PTS(:,1) +1i*PTS(:,2)
Interp_contour3  = interpolation1(my_contour3,scale,M,str)
histograms3 = extract_features(Interp_contour3,M,R,F,N2,s); % Histograms / Shape Context descriptors of 12x5 bins
    for i=1:M
    hist1=reshape(histograms3(:,:,i),1,[]);; % Reshaped single descriptor of 60D for each shape context
    H3=vertcat(H3,hist1); % 60D distribution for all points
    end
    figure,
    imagesc(~histograms3(:,:,i)');colormap gray


figure; imshow(X), title(' shapes for matching'),hold on,
plot(real(my_contour1),imag(my_contour1),'r.-');
plot(real(my_contour2)+x,imag(my_contour2)+x,'b.-');  
plot(real(my_contour3)+2*x,imag(my_contour3)+2*x,'g.-');  
drawnow

figure; imshow(X), title(' shapes after resampling'),hold on,
plot(real(Interp_contour1),imag(Interp_contour1),'r.-');
plot(real(Interp_contour2)+x,imag(Interp_contour2)+x,'b.-');
plot(real(Interp_contour3)+2*x,imag(Interp_contour3)+2*x,'g.-');
drawnow

% %% From the full body silhouette crop the "head+thorax" part
% % crop=imcrop(gcf)                       %# Crop the silhouette
% % imwrite(crop,'avatar.png')             %# Save the cropped silhouette
% % J=imread('avatar.png');                %# Load the silhouette data
% [PTS]=contour(BW)                       %# NN approach to get ordered shape contour
% my_contour1 = PTS(:,1) +1i*PTS(:,2)
% Interp_contour1  = interpolation1(my_contour1,scale,M,str)
% 
% %%       SHAPES FOR MATCHING 
% figure; imshow(X), title(' shapes for matching'),hold on,
% plot(real(my_contour1),imag(my_contour1),'r.-');
% 
% 
% %%      SHAPES AFTER RESAMPLING
% figure; imshow(X), title(' shapes after resampling'),hold on,
% plot(real(Interp_contour1),imag(Interp_contour1),'r.-');



