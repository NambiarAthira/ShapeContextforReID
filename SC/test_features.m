clear all, close all

scale = 10; M = 40; str = '';
R = 100; F = 12; N = 5; s = 'l';

% ----------------- Calculo de contornos interpolados ------------------ %

X = 255*ones(300,300);
figure; imshow(X),
[cp,lp,~,~] = get_points(X,[]);

my_contour = cp +1i*lp;

Interp_contour  = interpolation1(my_contour,scale,M,str);

%size(Interp_contour1) pode ser necessário se não o souberes graças à
%interpolação

%Código que escreve os pontos do contorno e respectivo índice dos pontos 
%perto deles (útil para perceber melhor os que estarão incluídos no
%círculo)
figure; imshow(X), hold on,
plot(real(Interp_contour),imag(Interp_contour),'r.');
label_points(Interp_contour, M);
hold off

%Colocar em cima de cada ponto de Interp_contour o disco de raio R,
%separado em F fatias e cada uma das fatias em N partes.
%Contar o número de pontos em cada uma das fatias obtendo um histograma 
%(matriz) FxN para cada ponto. A divisão espacial de cada uma das fatias
%pode ser uniforme ou logarítmica.
histograms = extract_features(Interp_contour,M,R,F,N,s);

%Função que desenha o "shape context" no ponto clicado pelo rato para 
%testar a correcta obtenção das features
draw_shapeContext(Interp_contour,M,R,F,N,s);