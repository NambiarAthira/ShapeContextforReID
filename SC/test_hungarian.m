%load hung_test.mat

n=400;
A=rand(n);

tic
[assignment, cost1] = munkres(A);
time1 = toc
tic
[assign, cost] = hungarian_alg(A);
time2 = toc
tic
[assig, cost2] = simplif_hung(A);
time3 = toc

m = assignment - assign;
c = cost1 - cost;
c2 = cost2 - cost;