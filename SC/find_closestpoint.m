%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Encontra o índice do array de pontos ("contourp") que identifica o ponto 
% mais próximo do ponto P = (x, y).
% contourp é o vector de pontos 
% np é o número de pontos
% x é a coordenada dos xx do ponto a procurar
% y é a coordenada dos yy do ponto a procurar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ind = find_closestpoint(contourp, np, x, y)
dist1 = euclid_distance(real(contourp(1)), imag(contourp(1)), x, y);
ind = 1;
for k = 2:np
    dist2 = euclid_distance(real(contourp(k)), imag(contourp(k)), x, y);
    if dist2 < dist1
        ind = k;
        dist1 = dist2;
    end
end