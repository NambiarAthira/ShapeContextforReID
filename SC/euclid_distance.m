%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Devolve a distância euclidiana entre dois pontos: (x1, y1) e (x2, y2)
% x1 coordenada xx do primeiro ponto
% y1 coordenada yy do primeiro ponto
% x2 coordenada xx do segundo ponto
% y2 coordenada yy do segundo ponto
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dist = euclid_distance(x1, y1, x2, y2)
dist = sqrt((x1 - x2)^2 + (y1 - y2)^2);