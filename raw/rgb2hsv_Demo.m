%function rgb2hsv_Demo()
clc;
close all;
clear all;

% Preallocate arrays to hold the means of all the images.
hImage_Mean = zeros(1, 1);
sImage_Mean = zeros(1, 1);
vImage_Mean = zeros(1, 1);

rgbImage= imread ('C:\HSV_MoG\MoG\DB4\3\3\rgbframe030.png')
figure,
subplot(3, 3, 1);
imshow(rgbImage);
caption = sprintf('Original Color Image\n%s');
title(caption, 'FontSize', 10)

% Convert to floating point so it does the calculations correctly.
% Also needs to be normalized to 0-1.
rgbFloating = double(rgbImage) / 255.0


% Compute hsv image
hsvImage = rgb2hsv(rgbImage);
% H image:
hImage = hsvImage(:,:,1);
subplot(3, 3, 4);
imshow(hImage, []); % Display the image.
% Compute mean
hImage_Mean= mean(hImage(:));
caption = sprintf('Hue Image. Mean = %6.2f', hImage_Mean);
title(caption, 'FontSize', 10);
% Compute and display the histogram for the H image.
% histogramDouble(hImage, 7, 'Histogram of Hue Image');
        minValue = min(min(hImage));
        maxValue = max(max(hImage));
        range = maxValue - minValue;
        dblImage = (hImage - minValue) / range;
        % Check to verify that range is now 0-1.
        % minValueNorm = min(min(dblImage));
        % maxValueNorm = max(max(dblImage));

        % Let's get its histogram into 256 bins.
        [pixelCount grayLevels] = imhist(dblImage, 256);

        % Let's suppress the zero bin because it's always so high.
        pixelCount(1) = 0;

        % But now grayLevelsD goes from 0 to 1.
        % We want it to go from the original range, so we need to scale.
        originalDoubleGrayLevels = range * grayLevels + minValue;

        subplot(3, 3, 5);
        plot(originalDoubleGrayLevels, pixelCount);
        title(caption, 'FontSize', 10);
        % Scale x axis manually.
        xlim([originalDoubleGrayLevels(1) originalDoubleGrayLevels(end)])



% S image:
sImage = hsvImage(:,:,2);
subplot(3,3,6);
imshow(sImage); % Display the image.
% Compute mean
sImage_Mean= mean(sImage(:));
caption = sprintf('Value Image. Mean = %6.2f', vImage_Mean);
title(caption, 'FontSize', 10);
% Compute and display the histogram for the H image.
% histogramDouble(hImage, 7, 'Histogram of Hue Image');
        minValue = min(min(sImage));
        maxValue = max(max(sImage));
        range = maxValue - minValue;
        dblImage = (sImage - minValue) / range;
        % Check to verify that range is now 0-1.
        % minValueNorm = min(min(dblImage));
        % maxValueNorm = max(max(dblImage));

        % Let's get its histogram into 256 bins.
        [pixelCount grayLevels] = imhist(dblImage, 256);

        % Let's suppress the zero bin because it's always so high.
        pixelCount(1) = 0;

        % But now grayLevelsD goes from 0 to 1.
        % We want it to go from the original range, so we need to scale.
        originalDoubleGrayLevels = range * grayLevels + minValue;

        subplot(3, 3, 7);
        plot(originalDoubleGrayLevels, pixelCount);
        title(caption, 'FontSize', 10);
        % Scale x axis manually.
        xlim([originalDoubleGrayLevels(1) originalDoubleGrayLevels(end)])



% V image:
vImage = hsvImage(:,:,3);
subplot(3, 3, 8);
imshow(vImage,[]); % Display the image.
% numberOfImagesProcessed = numberOfImagesProcessed + 1;
% Compute mean
vImage_Mean = mean(vImage(:));
caption = sprintf('Value Image. Mean = %6.2f', vImage_Mean);
title(caption, 'FontSize', 10);
% Compute and display the histogram for the V image.
% Compute and display the histogram for the H image.
% histogramDouble(vImage, 7, 'Histogram of Hue Image');
        minValue = min(min(vImage));
        maxValue = max(max(vImage));
        range = maxValue - minValue;
        dblImage = (vImage - minValue) / range;
        % Check to verify that range is now 0-1.
        % minValueNorm = min(min(dblImage));
        % maxValueNorm = max(max(dblImage));

        % Let's get its histogram into 256 bins.
        [pixelCount grayLevels] = imhist(dblImage, 256);

        % Let's suppress the zero bin because it's always so high.
        pixelCount(1) = 0;

        % But now grayLevelsD goes from 0 to 1.
        % We want it to go from the original range, so we need to scale.
        originalDoubleGrayLevels = range * grayLevels + minValue;

        subplot(3, 3, 9);
        plot(originalDoubleGrayLevels, pixelCount);
        title(caption, 'FontSize', 10);
        % Scale x axis manually.
        xlim([originalDoubleGrayLevels(1) originalDoubleGrayLevels(end)])


