clc
close all;
clear all;
nClus = 10;             % No of clusters(K)
I1=imread('SC1_1.png'); %  Person1 silhouette
[X1] = histo(I1)
size(X1)
X=X1'
size(X1)
size(X)
[pc,score] = princomp(X)
size(score)
plot(score(:,1),score(:,2),'.'),  
% hold on;
% 
% I2=imread('SC5_2.png'); %  Person1 silhouette
% [X1] = histo(I2)
% size(X1)
% X=X1'
% size(X1)
% size(X)
% [pc,score] = princomp(X)
% size(score)
% plot(score(:,1),score(:,2),'.') 
% 
