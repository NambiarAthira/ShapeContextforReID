%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Encontra o �ndice do array de pontos ("contourp") que identifica o ponto 
% mais pr�ximo do ponto P = (x, y).
% contourp � o vector de pontos 
% np � o n�mero de pontos
% x � a coordenada dos xx do ponto a procurar
% y � a coordenada dos yy do ponto a procurar
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