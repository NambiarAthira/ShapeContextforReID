
clear all, close all

scale = 10; M = 40; str = '';
R = 100; F = 12; N2 = 5; s = 'l';

% ----------------- Calculo de contornos interpolados ------------------ %

X = 255*ones(300,300);
figure; imshow(X),
title('2 shapes for matching')
[cp,lp,N,but] = get_points(X,[]);

my_contour1 = cp +1i*lp;

Interp_contour1  = interpolation1(my_contour1,scale,M,str);

%size(Interp_contour1) pode ser necess�rio se n�o o souberes gra�as �
%interpola��o

[cp,lp,N,but] = get_points(X,[]);

my_contour2 = cp +1i*lp;

Interp_contour2  = interpolation1(my_contour2,scale,M,str);

[cp,lp,N,but] = get_points(X,[]);
my_contour3 = cp +1i*lp;

Interp_contour3  = interpolation1(my_contour3,scale,M,str);

figure; imshow(X), hold on,
plot(real(Interp_contour1),imag(Interp_contour1),'r.');
plot(real(Interp_contour2),imag(Interp_contour2),'b.');
plot(real(Interp_contour3),imag(Interp_contour3),'k.');drawnow
title('Corresponding points are matched!')
hold off
%Colocar em cima de cada ponto de Interp_contour o disco de raio R,
%separado em F fatias e cada uma das fatias em N partes.
%Contar o n�mero de pontos em cada uma das fatias obtendo um histograma 
%(matriz) FxN para cada ponto. A divis�o espacial de cada uma das fatias
%pode ser uniforme ou logar�tmica.
histograms1 = extract_features(Interp_contour1,M,R,F,N2,s);
histograms2 = extract_features(Interp_contour2,M,R,F,N2,s);
histograms3 = extract_features(Interp_contour3,M,R,F,N2,s);

% ------------------- end shape context ------------------------ %
% 
% %Dados dois conjuntos de histogramas calcular a matriz de custo
% cm = cost_matrix(M,histograms1,histograms2);
% 
% %Dada a matriz de custo aplicar o algoritmo h�ngaro para resolver a
% %correspond�ncia entre os dois conjuntos de pontos.
% %Threshold???
% %cm2 = elim_min(cm); %esta fun��o pode n�o eliminar todos os m�nimos
% %repetidos. Em alternativa o problema da diferen�a entre o munkres e o
% %simplif_hung pode ser outro
% %TESTE
% savefile = 'hung_test.mat';
% save(savefile, 'cm');
% 
% tic
% [assignment, cost1] = munkres(cm);
% time1=toc
% %[assign, cost] = simplif_hung(cm);
% tic
% [assign, cost] = hungarian_alg(cm);
% time2=toc
% 
% m = assignment - assign;
% c = cost1 - cost;
% 
% %Dada a correspond�ncia obtida tra�ar linhas entre os pontos
% %correspondentes
% show_correspond(Interp_contour1, Interp_contour2, assignment, M);


figure; 
for i=1:40
subplot(4,10,i)
imagesc(~histograms1(:,:,i)');colormap gray
end
title('histograms of each point of shape1','FontWeight','bold')


figure;
for i=1:40
subplot(4,10,i)
imagesc(~histograms2(:,:,i)');colormap('gray')
end
title('histograms of each point of shape2','FontWeight','bold')


figure;
for i=1:40
subplot(4,10,i)
imagesc(~histograms3(:,:,i)');colormap('gray')
end
title('histograms of each point of shape3','FontWeight','bold')
