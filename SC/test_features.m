clear all, close all

scale = 10; M = 40; str = '';
R = 100; F = 12; N = 5; s = 'l';

% ----------------- Calculo de contornos interpolados ------------------ %

X = 255*ones(300,300);
figure; imshow(X),
[cp,lp,~,~] = get_points(X,[]);

my_contour = cp +1i*lp;

Interp_contour  = interpolation1(my_contour,scale,M,str);

%size(Interp_contour1) pode ser necess�rio se n�o o souberes gra�as �
%interpola��o

%C�digo que escreve os pontos do contorno e respectivo �ndice dos pontos 
%perto deles (�til para perceber melhor os que estar�o inclu�dos no
%c�rculo)
figure; imshow(X), hold on,
plot(real(Interp_contour),imag(Interp_contour),'r.');
label_points(Interp_contour, M);
hold off

%Colocar em cima de cada ponto de Interp_contour o disco de raio R,
%separado em F fatias e cada uma das fatias em N partes.
%Contar o n�mero de pontos em cada uma das fatias obtendo um histograma 
%(matriz) FxN para cada ponto. A divis�o espacial de cada uma das fatias
%pode ser uniforme ou logar�tmica.
histograms = extract_features(Interp_contour,M,R,F,N,s);

%Fun��o que desenha o "shape context" no ponto clicado pelo rato para 
%testar a correcta obten��o das features
draw_shapeContext(Interp_contour,M,R,F,N,s);