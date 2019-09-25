%% Regression of Shape Context with Biometric Features:

%% SC using the codes of Jacinto (functions of SC,COST,HUNGARIAN). 
% And does it for a series of automatically cropped top silhouettes (by means of selecting the minima point in neck,
% and thus get the chest point to crop).
% Then, the normalization of the height are carried out! (with a tolerance in BB)

clc;
close all;
clear all;

%%cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%% REAL SILHOUETTES (test)
%% Initialization
scale = 10;
M =40;   % the number of the contour points
str = '';
F = 12;   % Angular bins
N2 = 5;   % Radial bins
s = 'l';  % the type of division in "f": 'u' (uniforme) ou 'l' (logarítmica)
PTS=[];
S = struct(); 
S1=struct();
S2=struct(); 
addpath(genpath('C:\Users\AthiraNambiar\Desktop\PHD\codes\2014\ShapeContext\SC\SC'));

%% Get the avatar silhouette
 cd ('C:\Users\AthiraNambiar\Desktop\PHD\codes\2014\ShapeContext\SC\SC\examples_Human silhouettes\av')
%  cd ('C:\Users\AthiraNambiar\Desktop\PHD\codes\2014\ShapeContext\SC\SC\examples_Human silhouettes\avatars_togo')
imagefiles = dir('*.png');                  % get list of rgb .png files in this directory
noOfFrames = length(imagefiles)             % The structure of images details
images = strtrim(cellstr( num2str((1:noOfFrames)','person%d') )); 

for i=1:noOfFrames 
        im=imagefiles(i).name
        I=imread(im);   % Load the image data

 %% Get the Shape context histogram
        level = graythresh(I);                      % Compute an appropriate threshold
        BW = im2bw(I,level);                        % Convert grayscale to binary
        len=[];
        for idx=1:size(BW,1)
            len(idx)=size(find(BW(idx,:)==0),2);
        end
%        plot(1:size(BW,1),len)
        [mini]=filterr(len);
        [r,c]=min(abs(mini-round(size(BW,1)*7/12)));
        mini=mini(c);
        
%         figure, 
%         subplot(1,3,1),imshow(BW)
        crop_pt_bottom=round(mini+(mini- find(len,1))/2);
        crop_pt_top=round(find(len,1));
        vert_bound=find(imcomplement(BW(crop_pt_bottom,:)));
        crop_left=vert_bound(1);
        crop_right=vert_bound(end);
        width=crop_right-crop_left;
        height=crop_pt_bottom-crop_pt_top;
        tole=round(height*10/100);
%         I2=imcrop(BW,[crop_left,crop_pt_top,width,height]);
        I2=imcrop(BW,[crop_left-2*tole,crop_pt_top-tole,width+4*tole,height+2*tole]);
%         subplot(1,3,2),imshow(I2)
        sizI=size(I2);
        I = imresize(I2, [100 100*sizI(2)/sizI(1)]); % Normalize the height of the person 
%                                                      thus keeping the aspect ratio
%         subplot(1,3,3),
%          imshow(I)
%         title(['person:' num2str(im)])
%          im

        %% Get the Shape context histogram
        [PTS]=contour(I);                          % NN approach to get ordered shape contour
         my_contour1 = PTS(:,1) +1i*PTS(:,2);
        [Interp_contour1,perimeter] = interpolation1(my_contour1,scale,M,str);
        Xk=[real(Interp_contour1) imag(Interp_contour1)];
        R=perimeter/2;
%         R=100;
%         nsamp=size(Xk,1);
%         out_vec=zeros(1,nsamp);
%         r_inner=1/8;
%         r_outer=2;
%         nbins_theta=12;
%         nbins_r=5;
%         mean_dist_global=[];
        histograms = extract_features(Interp_contour1,M,R,F,N2,s); % Histograms / Shape Context descriptors of 12x5 bins
%         [histograms,mean_dist]=sc_compute(Xk',zeros(1,nsamp),mean_dist_global,nbins_theta,nbins_r,r_inner,r_outer,out_vec);
        histograms1=reshape(histograms,60,M);
        histograms2=reshape(histograms1,60*M,1);
        histo=histograms2';% 1 row vector histogram corresponding to SC of a person
        DATA(i,:)=histo;
        
        S.(images{i}) = histo;
%    
end
%% 

% Biometrics for the training set
    X= DATA
    N = size(X,1)
    samp_per=4; 
    samp_bio=9;
    samp_avatar=samp_per*samp_bio;
    no_avatars=size(X,1)/samp_avatar;
    
    E=[];
    RMSE=[];
    R2=[];
    AdjR2=[];
    bio1=[1    1  2  3  1  1  1  1  1 ]
    bio2=[1    1  1  1  2  3  1  1  1 ]
    bio3=[1    1  1  1  1  1 .5  2  1 ]
    bio4=[1.25 1  1  1  1  1  1  1  1 ]
    bio5=[1    1  1  1  1  1  1  1  1.25 ]
    y1=reshape(repmat(bio1,samp_per,1),samp_avatar,1) % Size_of_bio=9
    Y1=repmat(y1,no_avatars,1)
    y2=reshape(repmat(bio2,samp_per,1),samp_avatar,1)
    Y2=repmat(y2,no_avatars,1)
    y3=reshape(repmat(bio3,samp_per,1),samp_avatar,1)
    Y3=repmat(y3,no_avatars,1)
    y4=reshape(repmat(bio4,samp_per,1),samp_avatar,1)
    Y4=repmat(y4,no_avatars,1)
    y5=reshape(repmat(bio5,samp_per,1),samp_avatar,1);
    Y5=repmat(y5,no_avatars,1)
    Y=[Y1 Y2 Y3 Y4 Y5];
    
    
%CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC    

%% REAL SILHOUETTES (test)


%% Initialization
% scale = 10;
% M = 40;   % the number of the contour points
% str = '';
% F = 12;   % Angular bins
% N2 = 5;   % Radial bins
% s = 'l';  % the type of division in "f": 'u' (uniforme) ou 'l' (logarítmica)
% PTS=[];
T = struct(); 
% T1=struct();
% T2=struct(); 
% addpath(genpath('C:\Users\AthiraNambiar\Desktop\PHD\codes\2014\ShapeContext\SC\SC'));

% X = 255*ones(500,500);
%% Get the avatar silhouette
 cd ('C:\Users\AthiraNambiar\Desktop\PHD\codes\2014\ShapeContext\SC\SC\examples_Human silhouettes\hu')
%  cd ('C:\Users\AthiraNambiar\Desktop\PHD\codes\2014\ShapeContext\SC\SC\examples_Human silhouettes\human_trial')
imagefiles2 = dir('*.jpg');                  % get list of rgb .jpg files in this directory
noOfFrames2 = length(imagefiles2)             % The structure of images details
images2 = strtrim(cellstr( num2str((1:noOfFrames2)','person%d') )); 

for i=1:noOfFrames2 
        %% Normalize the height of the person
%         close all
        im2=imagefiles2(i).name
        I2=imread(im2);   % Load the image data
%         sizI=size(I);
%         I = imresize(I, [100 100*sizI(2)/sizI(1)]); % Normalize the height of the person 
%                                                     % thus keeping the aspect ratio

 %% Get the Shape context histogram
        level2 = graythresh(I2);                      % Compute an appropriate threshold
        BW2 = im2bw(I2,level2);                        % Convert grayscale to binary
        originalBW2=~BW2 % Invert the binary image(to have a white background just the way did for AVATARS)
        se2 = strel('disk',9)
        BW2 = imclose(originalBW2,se2)
        imshow(BW2)
        len2=[];
        for idx2=1:size(BW2,1)
            len2(idx2)=size(find(BW2(idx2,:)==0),2);
        end
%        plot(1:size(BW,1),len)
        [mini2]=filterr(len2);
        [r2,c2]=min(abs(mini2-round(size(BW2,1)*7/12)));
        mini2=mini2(c2);
         if mini2> 2/3*(size(len2,2))
            mini2=1/2*(size(len2,2))
        end
        
%         figure, 
%         subplot(1,3,1),imshow(BW)
        crop_pt_bottom2=round(mini2+(mini2- find(len2,1))/2);
        crop_pt_top2=round(find(len2,1));
        vert_bound2=find(imcomplement(BW2(crop_pt_bottom2,:)));
        crop_left2=vert_bound2(1);
        crop_right2=vert_bound2(end);
        width2=crop_right2-crop_left2;
        height2=crop_pt_bottom2-crop_pt_top2;
        tole2=round(height2*10/100);
%         I2=imcrop(BW,[crop_left,crop_pt_top,width,height]);
         I22=imcrop(BW2,[crop_left2-5*tole2,crop_pt_top2-tole2,width2+8*tole2,height2+7*tole2]);
%         subplot(1,3,2),imshow(I2)
        sizI2=size(I22);
        I2 = imresize(I22, [100 100*sizI2(2)/sizI2(1)]); % Normalize the height of the person 
%                                                      thus keeping the aspect ratio
%         subplot(1,3,3),
%         imshow(I)
%         title(['person:' num2str(im)])
%          im

        %% Get the Shape context histogram
        [PTS2]=contour(I2);                          % NN approach to get ordered shape contour
         my_contour12 = PTS2(:,1) +1i*PTS2(:,2);
        [Interp_contour12,perimeter2] = interpolation1(my_contour12,scale,M,str);
        Xk2=[real(Interp_contour12) imag(Interp_contour12)];
        R2=perimeter/2;
%         R=100;
%         nsamp=size(Xk,1);
%         out_vec=zeros(1,nsamp);
%         r_inner=1/8;
%         r_outer=2;
%         nbins_theta=12;
%         nbins_r=5;
%         mean_dist_global=[];
         histograms2 = extract_features(Interp_contour12,M,R,F,N2,s); % Histograms / Shape Context descriptors of 12x5 bins
%         [histograms,mean_dist]=sc_compute(Xk',zeros(1,nsamp),mean_dist_global,nbins_theta,nbins_r,r_inner,r_outer,out_vec);
        histograms12=reshape(histograms2,60,M);
        histograms22=reshape(histograms12,60*M,1);
        histo2=histograms22';% 1 row vector histogram corresponding to SC of a person
        DATA2(i,:)=histo2;
        
        T.(images2{i}) = histo2;
               
    end

   % Biometrics for the training set
    A= DATA2
  
  
    
    
    
 %% Regression Model

 y=[];  
 
   for k=1:5 % No of biometrics
    %     VARI=[];
           
     for i = 1:N
               
            %% PRINCIPLE COMPONENT REGRESSION
            x_train=X
            y_train=Y(:,k)
            x_test=A
%             y_test=B(:,k)

            [PCALoadings,PCAScores,PCAVar] = pca(x_train)
            betaPCR = regress(y_train, PCAScores(:,1:60))
            betaPCR = PCALoadings(:,1:60)*betaPCR;
            betaPCR = [mean(y_train) - mean(x_train)*betaPCR; betaPCR];
            yfitPCR = [ones(1,1) x_test]*betaPCR; %% Predicted BIOMETRICS!!!
%             yfitPCR = ([x_test]*betaPCR);
%             y_test=abs(yfitPCR)  %% Predicted BIOMETRICS!!!

     end
   y=[y y_test]
   end
  
   
   
%% Retrieval of people   
% load('boxplot_affine_10.mat', 'y')
boxplot(y)
XT = [1 2 3 4 5];
xlab={'neckness'  'chestSize'  'BodySize'  'HeadWidth'  'HeadLength'}
set(gca,'XTickLabel',xlab,'XTick',XT)
xlabel('Biometrics predicted in the real world data')
ylabel('Distribution of biometric data')

   

%% cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%%% Nearest Neighbour Test to verify Person Re-identification: (FIRST RANK)
% cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
% data=y;
% N=size(data,1);
% output=[];
% success=0;
% failure=0;
% 
% 
% for i = 1:N
%     test=data(i,:)
%     selection=zeros(N,1) 
%     selection(i)=1
%     train=data(find(selection==0),:)
%         
%     label=[repmat(1:length(selection)/4,4,1)];
%     group=reshape(label,[],1)
%     trainLabel=group(find(selection==0))
%     testLabel=group(find(selection==1))
%     train_group=testLabel
%        
%     test_group = knnclassify(test,train,trainLabel,1);
%     output= vertcat(output,[train_group test_group])
%     
%     if(train_group == test_group)
%         success=success+1
%     else
%         failure=failure+1
%     end
% end
%   accuracy= (success/(success+failure))*100
%   
  

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
% % %% Cumulative Rank  (CMC CURVE)
% 
% data=y
% cumulative_rank=zeros(length(y)/4,1)
% N=length(y)
% 
% for i=1:N
%     test=data(i,:)
%     selection=zeros(N,1) 
%     selection(i)=1
%     train=data(find(selection==0),:)
%         
%     label=[repmat(1:80,4,1)];
%     group=reshape(label,[],1)
%     trainLabel=group(find(selection==0))
%     testLabel=group(find(selection==1))
%         
%     scores = zeros(size(train,1),1)
%     for k=1:size(train,1)
%       
%              scores(k) = euclidean(train(k,:),test)
%      end
%            
%    % rcv_class = knnclassify(test,train,train_group,1)
%     [val,ind] = sort(scores)
%     
%      A=trainLabel(ind);
%     [uA a b] = unique(A,'first');
%     C= A(sort(a));
%         for q=1:10
%             if (C(q)==testLabel)
%                     cumulative_rank(q) = cumulative_rank(q) + 1;
%                     break;
%             end
%         end
%     
% %      for q=1:length(ind)
% %        
% %         if (trainLabel(ind(q)==testLabel))
% %         
% %             cumulative_rank(q) = cumulative_rank(q) + 1
% %             break;
% %             
% %         end
% %         
% %     end
%          
%     
% end
% cumulative_rank
% 
% 
% 
% cum_rank = zeros(length(cumulative_rank),1);
% 
% for i=1:length(cumulative_rank)
%    
%     cum_rank(i) = sum(cumulative_rank(1:i))/size(cumulative_rank,1)
%     
% end
% 
% x = 1:length(cumulative_rank);
% cum_rank = (cum_rank*100)/max(cum_rank)
% 
% figure,
% handle = plot(x,cum_rank,'r*-');
% xlabel('Cumulative Rank');
% ylabel('Recognition Rate (%)');
% xlim([1 20]) 
% ylim([0 100]) 
% 
% grid
% 
% 

