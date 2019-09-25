%% Cumulative Rank
cd('C:\Users\AthiraNambiar\Desktop\PHD\codes\2014\ShapeContext\SC\SC\MY_TRIALS\BipartiteMatchings\BIPARTITE_CROP_RESULT')
load('output_class.mat', 'y')
load('output_label.mat', 'label')
load('output_costShapes.mat', 'cost_Shapes')


no_people=54  % No of people -- [54 or 20]
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




