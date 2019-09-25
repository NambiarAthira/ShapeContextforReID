%% REAL WORLD DATA

%% SC using the codes of Jacinto (functions of SC,COST,HUNGARIAN). 
% And does it for a series of automatically cropped top silhouettes (by means of selecting the minima point in neck,
% and thus get the chest point to crop).
% Then, the normalization of the height are carried out! (with a  tolerance in BB)

clc;
close all;
clear all;

%% Initialization
scale = 10;
M = 40;   % the number of the contour points
str = '';
F = 12;   % Angular bins
N2 = 5;   % Radial bins
s = 'l';  % the type of division in "f": 'u' (uniforme) ou 'l' (logarítmica)
PTS=[];
S = struct(); 
S1=struct();
S2=struct(); 
addpath(genpath('C:\Users\AthiraNambiar\Desktop\PHD\codes\2014\ShapeContext\SC\SC'));

% X = 255*ones(500,500);
%% Get the avatar silhouette
 cd ('C:\Users\AthiraNambiar\Desktop\PHD\codes\2014\ShapeContext\SC\SC\examples_Human silhouettes\hu')
% cd ('C:\Users\AthiraNambiar\Desktop\PHD\codes\2014\ShapeContext\SC\SC\examples_Human silhouettes\probs')
imagefiles = dir('*.jpg');                  % get list of rgb .jpg files in this directory
noOfFrames = length(imagefiles)             % The structure of images details
images = strtrim(cellstr( num2str((1:noOfFrames)','person%d') )); 

for i=1:noOfFrames 
        %% Normalize the height of the person
%         close all
        im=imagefiles(i).name
        I=imread(im);   % Load the image data
%         sizI=size(I);
%         I = imresize(I, [100 100*sizI(2)/sizI(1)]); % Normalize the height of the person 
%                                                     % thus keeping the aspect ratio

 %% Get the Shape context histogram
        level = graythresh(I);                      % Compute an appropriate threshold
        BW = im2bw(I,level);                        % Convert grayscale to binary
        originalBW=~BW % Invert the binary image(to have a white background just the way did for AVATARS)
        se = strel('disk',9)
        BW = imclose(originalBW,se)
%         imshow(BW)
        len=[];
        for idx=1:size(BW,1)
            len(idx)=size(find(BW(idx,:)==0),2);
        end
%        plot(1:size(BW,1),len)
        [mini]=filterr(len);
        [r,c]=min(abs(mini-round(size(BW,1)*7/12)));
        mini=mini(c);
         if mini> 2/3*(size(len,2))
            mini=1/2*(size(len,2))
        end
        
%         figure, 
%         subplot(1,3,1),
        crop_pt_bottom=round(mini+(mini- find(len,1))/2);
        crop_pt_top=round(find(len,1));
        vert_bound=find(imcomplement(BW(crop_pt_bottom,:)));
        crop_left=vert_bound(1);
        crop_right=vert_bound(end);
        width=crop_right-crop_left;
        height=crop_pt_bottom-crop_pt_top;
        tole=round(height*10/100);
%         imshow(BW)
%         I2=imcrop(BW,[crop_left,crop_pt_top,width,height]);
        I2=imcrop(BW,[crop_left-5*tole,crop_pt_top-tole,width+8*tole,height+7*tole]);
%         subplot(1,3,2),imshow(I2)
        sizI=size(I2);
        I = imresize(I2, [100 100*sizI(2)/sizI(1)]); % Normalize the height of the person 
%                                                      thus keeping the aspect ratio
%         subplot(1,3,3),
          imshow(I)
          sizi=size(I)
%         title(['person:' num2str(im)])
%          im

        %% Get the Shape context histogram
        [PTS]=contour(I);                          % NN approach to get ordered shape contour
         my_contour1 = PTS(:,1) +1i*PTS(:,2);
        [Interp_contour1,perimeter] = interpolation1(my_contour1,scale,M,str);
        Xk=[real(Interp_contour1) imag(Interp_contour1)];
        R=perimeter/2;
%         
         histograms = extract_features(Interp_contour1,M,R,F,N2,s); % Histograms / Shape Context descriptors of 12x5 bins
%         [histograms,mean_dist]=sc_compute(Xk',zeros(1,nsamp),mean_dist_global,nbins_theta,nbins_r,r_inner,r_outer,out_vec);
        
            %         histograms1=reshape(histograms,60,M);
            %         histograms2=reshape(histograms1,60*M,1);
            %         histo=histograms2';% 1 row vector histogram corresponding to SC of a person
            %         DATA(i,:)=histo;
            %         S.(images{i}) = histo;
               
            S.(images{i}) = histograms;
            if mod(i,2)==1
                S1.(images{i})= histograms;         % Trainset
            else
                S2.(images{i})= histograms;         % Testset
            end
               
    end

     
        
 %% Cost matching of SC 
noPersons=noOfFrames/4;     % Number of persons
x=2*noPersons;              % Number of frames in each test/train set per person
cost_Shapes=[];


for i=1:2:noOfFrames
    for j=2:2:noOfFrames
       histogram1=S1.(images{i});                  % histogram of first shape
       histogram2=S2.(images{j});                  % histogram of second shape
           cm = cost_matrix(M,histogram1,histogram2);  % Returns an array of cost that compares two sets of two different SC.
%          cm=hist_cost_2(histogram1,histogram2);
%           [assign, cost] = hungarian_alg(cm);
        [assign, cost]=hungarian(cm);              % Malik code on Hungarian algorithm           
       cost_Shapes=[cost_Shapes;cost];
    end
end
cost_Shapes=reshape(cost_Shapes,x,x)

label=[repmat(1:noPersons,2,1)];
label=reshape(label,[],1);
[vals, kindex] = min(cost_Shapes);
y = label(kindex)

%% Confusion Matrix
[C,order] = confusionmat(label,y)
CM=diag(1./sum(C,2))*C*100
%Visualize the Confusion Matrix
Accuracy= trace(C)/length(label)
fprintf('NN Classifier Accuracy: %.2f%%\n', 100*Accuracy)
imagesc(CM);
colorbar; 
title(['CM with', num2str(noPersons) 'people, NN Classifier Accuracy:', num2str(100*Accuracy)])
noPersons       

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
%% Cumulative Rank

no_people=noPersons  % No of people
no_samples=2 % No of samples in test per person
smples_test=no_people*no_samples
[val,ind] = sort(cost_Shapes) % Sort the cost_Shapes
C=[]         % Initialization of the matrix C
A=ceil(ind/no_samples) % ind replaced by corresponding labels of identified person
cumulative_rank = zeros(no_people,1) % Initialize the cumulative_rank vector

% Get the compressed matrix with unique id in the training.(Each row corresponding to a rank)
for i=1:smples_test
    c=unique(A(:,i),'stable')
    C=[C c] % Concatenate the best reidentification for each column as a matrix.
end
C

% Get the re-identification rank score for each trial
for z=1:size(C,1)% Rank no
    for i=1:size(C,2) % no of test samples
        if (C(z,i)==label(i))
            cumulative_rank(z) = cumulative_rank(z) + 1
%             break;
        end
    end
     
end
cumulative_rank

cum_rank = zeros(size(cumulative_rank))

for i=1:length(cumulative_rank)
       cum_rank(i) = sum(cumulative_rank(1:i))
end
cmc=cum_rank/max(cum_rank)*100

x = 1:no_people
y= cmc


figure,
handle = plot(x,y,'r*-');
xlabel('Cumulative Rank');
ylabel('Recognition Rate (%)');
xlim([1 20])
ylim([50 100]) 

grid



