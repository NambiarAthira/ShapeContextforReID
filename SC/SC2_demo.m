%% Demo of SC2.m .Here, applied on person silhouettes and get their Shape Context 

%%


clc,
close all,
clear all, 

scale = 10; M = 40; str = '';x=70;
R = 100; F = 12; N2 = 5; s = 'l';
H1=[];H2=[];H3=[];H4=[];
PTS=[];

X = 255*ones(300,300);
% figure; imshow(X), hold on;
% ----------------- Calculo de contornos interpolados ------------------ %

% PERSON1
I=imread('SC1_1.png')
[PTS]=contour(I)
my_contour1 = PTS(:,1) +1i*PTS(:,2)
Interp_contour1  = interpolation1(my_contour1,scale,M,str);


% PERSON2
I2 = imread('SC2_1.png')
[PTS2]=contour(I2)
my_contour2 = PTS2(:,1) +1i*PTS2(:,2)
Interp_contour2  = interpolation1(my_contour2,scale,M,str);

% % PERSON3
% I3 = imread('SC4_2.png')
% [PTS3]=contour(I3)
% my_contour3 = PTS3(:,1) +1i*PTS3(:,2)
% Interp_contour3  = interpolation1(my_contour3,scale,M,str);
% 
% % PERSON4
% I4 = imread('SC5_2.png')
% [PTS4]=contour(I4)
% my_contour4 = PTS4(:,1) +1i*PTS4(:,2)
% Interp_contour4  = interpolation1(my_contour4,scale,M,str);

figure; imshow(X), title(' shapes for matching'),hold on,
plot(real(my_contour1),imag(my_contour1),'r.-');
plot(real(my_contour2)+x,imag(my_contour2)+x,'b.-');  
% plot(real(my_contour3)+2*x,imag(my_contour3)+2*x,'k.-');  
% plot(real(my_contour4)+3*x,imag(my_contour4)+3*x,'r.-'); 
drawnow

figure; imshow(X), title(' shapes after resampling'),hold on,
plot(real(Interp_contour1),imag(Interp_contour1),'r.-');
plot(real(Interp_contour2)+x,imag(Interp_contour2)+x,'b.-');
% plot(real(Interp_contour3)+2*x,imag(Interp_contour3)+2*x,'k.-');
% plot(real(Interp_contour4)+3*x,imag(Interp_contour4)+3*x,'r.-');
drawnow

% title('Corresponding points are matched!')
hold off
%Colocar em cima de cada ponto de Interp_contour o disco de raio R,
%separado em F fatias e cada uma das fatias em N partes.
%Contar o número de pontos em cada uma das fatias obtendo um histograma 
%(matriz) FxN para cada ponto. A divisão espacial de cada uma das fatias
%pode ser uniforme ou logarítmica.
histograms1 = extract_features(Interp_contour1,M,R,F,N2,s);
for i=1:40
hist1=reshape(histograms1(:,:,i),1,[])
H1=vertcat(H1,hist1)
end
histograms2 = extract_features(Interp_contour2,M,R,F,N2,s);
for i=1:40
hist2=reshape(histograms2(:,:,i),1,[])
H2=vertcat(H2,hist2)
end
% histograms3 = extract_features(Interp_contour3,M,R,F,N2,s);
% for i=1:40
% hist3=reshape(histograms3(:,:,i),1,[])
% H3=vertcat(H3,hist3)
% end
% histograms4 = extract_features(Interp_contour4,M,R,F,N2,s);
% for i=1:40
% hist4=reshape(histograms4(:,:,i),1,[])
% H4=vertcat(H4,hist4)
% end


%% Difference of sc histograms at different locations in a silhouette
%(corresponding to the sample Shape context figure in our paper)

% H1
% figure; 
% subplot(1,3,1)
% sc=reshape(H1(25,:),5,12) % at a position of 25st sample point
% sc_new=1-sc
% imagesc(sc_new);
% colormap gray
% xlabel('\theta')
% ylabel('logr')
% 
% subplot(1,3,2)
% sc=reshape(H1(26,:),5,12) % at a position of 26st sample point
% sc_new=1-sc
% imagesc(sc_new);
% colormap gray
% xlabel('\theta')
% ylabel('logr')
% 
% subplot(1,3,3)
% sc=reshape(H1(6,:),5,12) % at a position of 6th sample point
% sc_new=1-sc
% imagesc(sc_new);
% colormap gray
% xlabel('\theta')
% ylabel('logr')


%% Silhouette histogram of person, by flattening and concatenating the sc histograms

figure; 
subplot(1,4,1)
h1=1-H1
imagesc(~H1);
colormap gray
subplot(1,4,2)
h2=1-H2
imagesc(~H2);
colormap gray

% subplot(1,4,3)
% imagesc(~H3);colormap gray
% subplot(1,4,4)
% imagesc(~H4);colormap gray

% ------------------- end shape context ------------------------ %
% % 
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
show_correspond(Interp_contour1, Interp_contour2+100, assignment, M);

%% TO SHOW EACH POINT SHAPE CONTEXT HISTOGRAM OF 12*5
% figure; 
% for i=1:40
% subplot(4,10,i)
% imagesc(~histograms1(:,:,i)');colormap gray
% end
% title('histograms of each point of shape1','FontWeight','bold')
% 
% 
% figure;
% for i=1:40
% subplot(4,10,i)
% imagesc(~histograms2(:,:,i)');colormap('gray')
% end
% title('histograms of each point of shape2','FontWeight','bold')
% 
% figure;
% for i=1:40
% subplot(4,10,i)
% imagesc(~histograms3(:,:,i)');colormap('gray')
% end
% title('histograms of each point of shape3','FontWeight','bold')
