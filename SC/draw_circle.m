%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desenha um círculo com centro(cx, cy) e raio radius.
% cx é a coordenada xx do centro
% cy é a coordenada yy do centro
% radius é o raio do círculo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_circle(cx, cy, radius)
t = 0:0.01:2*pi;
x = (radius*cos(t))+cx;
y = (radius*sin(t))+cy;
plot(x,y);