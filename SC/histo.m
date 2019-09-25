function [X] = histo(I)

scale = 10; M = 40; str = '';
% R = 100;
F = 12; N2 = 5; s = 'l';H1=[];

[PTS]=contour(I);    %  Contour points  
my_contour1 = PTS(:,1) +1i*PTS(:,2);
[Interp_contour1,perimeter]= interpolation1(my_contour1,scale,M,str); % Resampled contour points
R=perimeter/4;
histograms1 = extract_features(Interp_contour1,M,R,F,N2,s); % Histograms / Shape Context descriptors of 12x5 bins
    for i=1:M
    hist1=reshape(histograms1(:,:,i),1,[]);; % Reshaped single descriptor of 60D for each shape context
    H1=vertcat(H1,hist1); % 60D distribution for all points
    end
X=H1'; % D x N matrix representing feature vectors by columns where 
       % D is the number of dimensions and N is the number of vectors.(60x40)
end
 
