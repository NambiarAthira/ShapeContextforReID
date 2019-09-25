%% Code for Nearest Neighbour Algorithm using Kullback Leiber Divergence

% With 3 persons' HSV Histograms; 100 trainset and 100 testset per person
%Author: Athira
%Date  : 07/03/2014
 
%% Reset

clc; 
close all;
clear all;

%% Read the input data
noPersons=10
train_person=100
test_person=100

S = struct();                                                         % The structure of images details
images = strtrim(cellstr( num2str((1:noPersons)','im%d') ));                  %'# field1,field2,...

for i=1:noPersons
    filename=strcat(int2str(i),'.mat');
    S.(images{i}) = load(filename);         % input each HSVhist data under each image folder
end

for i=1:noPersons
xtrain(:,(i-1)*train_person+1:i*train_person) = S.(images{i}).HSV_person(1:train_person,:)' ;  % 100 36-dimesional training points from class 1
end
training=xtrain';  

label=[repmat(1:noPersons,train_person,1)];
group=reshape(label,[],1);

for i=1:noPersons
xtest(:,(i-1)*test_person+1:i*test_person) = S.(images{i}).HSV_person(100+1:100+test_person,:)' ;  % 100 36-dimesional training points from class 1
end

sample=xtest';
class = knnclassify(sample, training, group,1) % Knn classify
[C,order] = confusionmat(group,class)
%Visualize the Confusion Matrix
imagesc(C);
colorbar;

Accuracy= trace(C)/length(group)
fprintf('NN Classifier Accuracy: %.2f%%\n', 100*Accuracy)


%   
% 
% 
