uA = [0.1,0.1]; uB =[-0.1,-0.1];
sA = 1000; sB = 100;
%% Trainset
xTr = vertcat( bsxfun(@plus,uA,randn(sA,2)),bsxfun(@plus,uB,randn(sB,2)));
yTr = vertcat(ones(sA,1),ones(sB,1)*-1);
%% Testset
xTe = vertcat( bsxfun(@plus,uA,randn(sA,2)),bsxfun(@plus,uB,randn(sB,2)));
yTe =  vertcat(ones(sA,1),ones(sB,1)*-1);
%% regression coefficient
w = xTr\yTr;
errRateNoConstantNotCentered_1 = mean(sign(xTe*w)==yTe)

uX = mean(xTr,1);
wCentered = bsxfun(@minus,xTr,uX)\yTr;
errRateCenteredNoConstant_1 = mean(sign( bsxfun(@minus,xTe,uX)*wCentered)==yTe)

wConstant = [ones(size(xTr,1),1) xTr]\yTr;
errRateConstantTerm_1 = mean(sign([ones(size(xTe,1),1) xTe]*wConstant)==yTe)

%% Athira TRIAL

% Not Centered, No Affine
w = xTr\yTr; % Regression coefficient
yPred=xTe*w;
errRateNoConstantNotCentered_2 = mean(sign(yPred)==yTe)

% Centered , No Affine
uX = mean(xTr,1);
uY = mean(yTr,1);
wCentered = bsxfun(@minus,xTr,uX)\yTr;
errRateCenteredNoConstant_2 = mean(sign((bsxfun(@minus,xTe,uX)*wCentered)+uY)==yTe)

% Not Centered, Affine
wConstant = [ones(size(xTr,1),1) xTr]\yTr;
errRateConstantTerm_2 = mean(sign([ones(size(xTe,1),1) xTe]*wConstant)==yTe)