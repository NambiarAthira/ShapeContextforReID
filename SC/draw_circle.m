%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desenha um c�rculo com centro(cx, cy) e raio radius.
% cx � a coordenada xx do centro
% cy � a coordenada yy do centro
% radius � o raio do c�rculo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_circle(cx, cy, radius)
t = 0:0.01:2*pi;
x = (radius*cos(t))+cx;
y = (radius*sin(t))+cy;
plot(x,y);