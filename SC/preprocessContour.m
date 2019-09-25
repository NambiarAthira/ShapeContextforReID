% Preprocessing for obtaining the human curvature 
I = imread('SC1_1.png')
BW1 = edge(I,'prewitt'); 
imshow(BW1)
[y,x]=find(BW1==1)
ix=[x y]

% Find the extreme points on the shape at the bottomlevel in order to find
% the centroid
find(y==max(y))
extrema=ix(r,:)
centroid=mean(extrema)




my_contour1 = x +1i*y;
Interp_contour1  = interplation1(my_contour1,scale,M,str);