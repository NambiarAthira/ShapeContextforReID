function [PTS]=contour3(I)

PTS=[];

% I=logical( I/255);     % Logical
I = edge(I,'canny'); %Edge detection
% figure; imshow(I),hold on;
 
[lin, col] = find(I>=1);   
x = [col lin] ;      % Co ordintes of binary edge

%% Ordering/ Resampling of points

[r]=find((x(:,2))==max(x(:,2))); % get the extreme points where whe need to open the shape
idx=min(r); % Minimum among the extremas a.k.a first point
y=x(idx,:);% Starting point
endpt=x(end,:);% Ending point
% plot(y(1),y(2),'ro')
% plot(y(1),y(2),'r.-')
% plot(endpt(1),endpt(2),'c*')
x(idx,:)=[]; % Remove the already marked point from the searching points
PTS=[PTS;y] ;

while y(1)~= endpt(1)|| y(2)~= endpt(2)
    [PI,D]=knnsearch(x,y,'k',1,'distance','euclidean');% Nearest neighbour of point y
    y=x(PI,:); % Neighbour point
%     plot(y(1),y(2),'ro')
%     plot(y(1),y(2),'r.-')
    PTS=[PTS;y]; % Ordered/ resampled points 
    x(PI,:)=[];% Remove already marked point from the searching points
end
end  
