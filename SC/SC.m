%% The basic function to draw image and  get the SC.

%%
clear all, close all

scale = 10; M = 40; str = '';
R = 100; F = 12; N2 = 5; s = 'l';

% ----------------- Calculo de contornos interpolados ------------------ %

X = 255*ones(300,300);
figure; imshow(X),

clc
[cp,lp,N,but] = get_points(X,[]);
my_contour1 = cp +1i*lp;
Interp_contour1  = interpolation1(my_contour1,scale,M,str);

%size(Interp_contour1) pode ser necessário se não o souberes graças à
%interpolação

%[cp,lp,N,but] = get_points(X,[]);
%my_contour2 = cp +1i*lp;
%Interp_contour2  = interpolation1(my_contour2,scale,M,str);

figure; imshow(X), hold on,
plot(real(Interp_contour1),imag(Interp_contour1),'r.');
%plot(real(Interp_contour2),imag(Interp_contour2),'b.'); drawnow
hold off
%Colocar em cima de cada ponto de Interp_contour o disco de raio R,
%separado em F fatias e cada uma das fatias em N partes.
%Contar o número de pontos em cada uma das fatias obtendo um histograma 
%(matriz) FxN para cada ponto. A divisão espacial de cada uma das fatias
%pode ser uniforme ou logarítmica.
histograms1 = extract_features(Interp_contour1,M,R,F,N2,s);
%histograms2 = extract_features(Interp_contour2,M,R,F,N2,s);

pause

% ------------------- end shape context ------------------------ %

%Dados dois conjuntos de histogramas calcular a matriz de custo
cm = cost_matrix(M,histograms1,histograms2);

%Dada a matriz de custo aplicar o algoritmo húngaro para resolver a
%correspondência entre os dois conjuntos de pontos.
%Threshold???
%cm2 = elim_min(cm); %esta função pode não eliminar todos os mínimos
%repetidos. Em alternativa o problema da diferença entre o munkres e o
%simplif_hung pode ser outro
%TESTE
savefile = 'hung_test.mat';
save(savefile, 'cm');

tic
[assignment, cost1] = munkres(cm);
time1=toc
%[assign, cost] = simplif_hung(cm);
tic
[assign, cost] = hungarian_alg(cm);
time2=toc

m = assignment - assign;
c = cost1 - cost;

%Dada a correspondência obtida traçar linhas entre os pontos
%correspondentes
show_correspond(Interp_contour1, Interp_contour2, assignment, M);

figure;
for i=1:10
subplot(2,5,i)
imagesc(histograms1(:,:,i)');colormap gray
end
