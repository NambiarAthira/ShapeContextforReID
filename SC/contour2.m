function [PTS]=contour2(I)

PTS=[];
tole=3;

% I=logical( I/255);     % Logical
I = edge(I,'canny'); %Edge detection
% figure; imshow(I),hold on;
 
[lin, col] = find(I>=1);   
x = [col lin] ;      % Co ordintes of binary edge

%% Ordering/ Resampling of points

[r]=find((x(:,2))==max(x(:,2))); % get the extreme points where whe need to open the shape
idx=min(r); % Minimum among the extremas a.k.a first point
y=x(idx,:);% Starting point
start_pt=y;
endpt=x(end,:);% Ending point
end_pt=endpt;
width=end_pt(1)-start_pt(1);
height=start_pt(2)-min(x(:,2));% Height
I2=imcrop(I,[start_pt(1)-tole start_pt(2)-height-tole width+tole height]);
% imshow(I2)
% figure,
asp_ratio=width/height;
I3 = imresize(I2, [200 200*asp_ratio]);
% imshow(I3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[lin, col] = find(I3>=1);   
x = [col lin] ;      % Co ordintes of binary edge
%Ordering/ Resampling of points

[r]=find((x(:,2))==max(x(:,2))); % get the extreme points where whe need to open the shape
idx=min(r); % Minimum among the extremas a.k.a first point
y=x(idx,:);% Starting point
endpt=x(end,:);% Ending point

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

