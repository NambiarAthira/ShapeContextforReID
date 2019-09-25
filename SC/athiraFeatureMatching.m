%%Code for feature matching (CMC curve is plotted against all gallery samples eg. 30 samples for 3 persons , its not what we need.
% We need it to be plotted against the no of users/persons which is 3 )
%Create Trainset and Test set and compare using euclidean and Mahalanobis distances
% Author: Athira
%Date: 30-01-2013
% 
% clc;
% close all;
% clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EUCLIDEAN DISTANCE
%Create Trainset and Test set 
% 1) Generating the trainset

G=struct();
 
images = strtrim(cellstr( num2str((1:3)','im%d') ));                  %'# field1,field2,...
for i=1:numel(images)
    G.(images{i}).selection=S.(images{i}).data(SSS.(images{i}).trimImages(1:end,end),:);
end
permu=randperm(size(G.im1.selection,1));
Trainset1=G.im1.selection(permu(1:10),:) 
Testset1=G.im1.selection(permu(11:15),:)

permu=randperm(size(G.im2.selection,1));
Trainset2=G.im2.selection(permu(1:10),:)
Testset2=G.im2.selection(permu(11:15),:) 

permu=randperm(size(G.im3.selection,1));
Trainset3=G.im3.selection(permu(1:10),:)
Testset3=G.im3.selection(permu(11:15),:) 

Trainset=[Trainset1;Trainset2;Trainset3]
Testset=[Testset1;Testset2;Testset3]

testData=Testset
trainData=Trainset

% Save data for future use
matfile = ['db-data.mat'];
save(matfile, 'trainData', 'testData');

load('C:\Users\AthiraNambiar\Desktop\PHD\codes\MoG\db-data.mat')
load('C:\Users\AthiraNambiar\Desktop\PHD\codes\MoG\db-labels.mat')

%D_eucl= pdist2(trainData,testData,'euclidean') % this will calculate the
% euclidean distance in single line :-)

for i=1:size(testData,1) % no of test sequences
for j=1:size(trainData,1)  % no of train sequences
euc_distances(j,i) = euclidean(trainData(j,:),testData(i,:));
end
end
save('euc_distances_HSV.mat','euc_distances')
% %Euclidean Matching
% 
% for i=1:size(Testset,1)
% for j=1:size(Trainset,1)
% euc_distances(j,i) = euclidean(Trainset(j,:),Testset(i,:));
% end
% end
% % save('euc_distances_pca.mat','euc_distances')
% % for i=1:20
% % for j=1:4
% % euc_distances(j,i) = euclidean(trainP(j,:),testP(i,:));
% % end
% % end
% % save('euc_distances_pca.mat','euc_distances')
% 
% euc_distances
% [r,c] = min(euc_distances)
% c'
%position_pq = [1:size(euc_distances,2);ii]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
%% Cumulative Rank

train = trainData;
test = testData;

% runSvdPca('.')
% load('C:\HSV_MoG\MoG\db-svdpca.mat')

num_users = size(train,1);

cumulative_rank = zeros(num_users,1);

for z=1:size(test,1)

    rcv = test(z,:);

    scores = zeros(num_users,1);

    tic;

    for i=1:num_users
    
        scores(i) = euclidean(train(i,:),rcv);
        
    end
    
    toc;
    
    [val,ind] = sort(scores);
    
    for i=1:length(ind)
       
        if (strcmp(trainLabels{ind(i)},testLabels{z}))
        
            cumulative_rank(i) = cumulative_rank(i) + 1;
            break;
            
        end
        
    end
    
end

cum_rank = zeros(num_users,1);

for i=1:length(cumulative_rank)
   
    cum_rank(i) = sum(cumulative_rank(1:i))/size(test,1);
    
end

x = 1:num_users;
cum_rank = cum_rank(1:num_users)*100;

figure,
handle = plot(x,cum_rank,'r*-');
xlabel('Cumulative Rank');
ylabel('Recognition Rate (%)');
ylim([0 100]) 

grid