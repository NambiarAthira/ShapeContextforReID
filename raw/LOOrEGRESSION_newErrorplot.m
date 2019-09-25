%% Regression of Shape Context with Biometric Features:

%% SC using the codes of Jacinto (functions of SC,COST,HUNGARIAN). 
% And does it for a series of automatically cropped top silhouettes (by means of selecting the minima point in neck,
% and thus get the chest point to crop).
% Then, the normalization of the height are carried out! (with a tolerance in BB)

clc;
close all;
clear all;

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


for k=1:5 % No of biometrics
    MU=[];
    VARI=[];
        
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
    
    bio1=[1  1  2  3  1  1  1  1  1 ]
    bio2=[1  1  1  1  2  3  1  1  1 ]
    bio3=[1  1  1  1  1  1  .5 2  1 ]
    bio4=[1.25 1 1 1  1  1  1  1  1 ]
    bio5=[1  1  1  1  1  1  1  1  1.25 ]
    y1=reshape(repmat(bio1,samp_per,1),samp_avatar,1) % Size_of_bio=9
    Y1=repmat(y1,no_avatars,1)
    y2=reshape(repmat(bio2,samp_per,1),samp_avatar,1)
    Y2=repmat(y2,no_avatars,1)
    y3=reshape(repmat(bio3,samp_per,1),samp_avatar,1)
    Y3=repmat(y3,no_avatars,1)
    y4=reshape(repmat(bio4,samp_per,1),samp_avatar,1)
    Y4=repmat(y4,no_avatars,1)
    y5=reshape(repmat(bio5,samp_per,1),samp_avatar,1)
    Y5=repmat(y5,no_avatars,1)
    Y=[Y1 Y2 Y3 Y4 Y5];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Regression Model

            for i = 1:N
%                 [train,test] = crossvalind('LeaveMOut',N,1);
%                 [r,c]=find(train==0);
             test= zeros(N,1)
             test(i)=1
             train= ones(N,1)
             train(i)=0


            %% PRINCIPLE COMPONENT REGRESSION

            x_train=X(find(train==1),:)
            y_train=Y(find(train==1),k)
            x_test=X(find(test==1),:)
            y_test=Y(find(test==1),k)

            [PCALoadings,PCAScores,PCAVar] = pca(x_train)
            betaPCR = regress(y_train, PCAScores(:,1:end))
            betaPCR = PCALoadings(:,1:end)*betaPCR;
            betaPCR = [mean(y_train) - mean(x_train)*betaPCR; betaPCR];
            yfitPCR = [ones(1,1) x_test]*betaPCR; %Predicted output
            e=abs((y_test-yfitPCR)) % Difference between true and predicted output
            E=[E;e];
            

            varpY=mean(((y_train-mean(y_train)).^2))% Population Varianceof y/ Average squared deviation of y from its mean
            varY=var(y_train)                 % Sample variance of y
            D=size(x_train,1)-2               % Degrees of Freedom
            rmse=sqrt(sum((y_test-yfitPCR).^2)/D);
            RSquare=1-((mean((y_test-yfitPCR).^2)/varpY))
            AdjRSquare=1-(rmse.^2/varY)
            
            RMSE=[RMSE;rmse];
            R2=[R2;RSquare];
            AdjR2=[AdjR2;AdjRSquare];

            end
%             [E RMSE R2 AdjR2]
E
            
    % NORMALIZE  THE PREDICTED VALUES CORRESPONDING TO THE RANGE OF THE
    % DATA USED FOR EVALUATION
    range=abs(max(y_train)-min(y_train))
    EE=E/range
    betaPCR=[];
    mu=mean(EE)
    MU=vertcat(MU,mu);
    vari=var(EE)
    VARI=vertcat(VARI,vari)
 

    MEAN(k)=mean(MU)
    VARIANCE(k)=mean(VARI)
    
    ROOTMEANSQERR(k)=mean(RMSE)
    R_SQUARE(k)=mean(R2)
    ADJ_R_SQUARE(k)=mean(AdjR2)

end

errorbar(MEAN,VARIANCE,'rx')
title('Errorplot for the regression Analysis')
XT = [1 2 3 4 5];
xlab={'neckness'  'chestSize'  'BodySize'  'HeadWidth'  'HeadLength'}
set(gca,'XTickLabel',xlab,'XTick',XT)
xlabel('Biometric') % x-axis label
ylabel('|TrueY-ObservedY|') % y-axis label













