%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Retorna uma matriz de custo que compara dois dois diferentes conjuntos de
% shape contexts
% np é o número de pontos a que foram aplicados os shape contexts - vamos 
% assumir que as duas figuras a comparar têm o mesmo número de pontos, o
% que faz sentido pois é feita uma interpolação para cada figura. Caso não
% se pudesse assumir o mesmo número de pontos era necessário acrescentar
% componentes com infinito (?) para a figura com menos pontos.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Returns an array of cost that compares two sets of two different 
% Shape contexts
% Np is the number of points that the shape contexts were applied - we 
% Assume comparing the two figures have the same number of points, 
% Which makes sense as an interpolation is made for each figure. if not 
% If you could take the same number of points was necessary to add 
% Components with infinity (?) For the figure with less points.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





function costm = cost_matrix(np, hist1, hist2)
costm = zeros(np, np);
for k = 1:np
    for l = 1:np
        costm(k, l) = costh(hist1(:,:,k), hist2(:,:,l)); %dados dois 
        %histogramas dos pontos k e l de diferentes figuras 
        %esta função devolve um valor que representa o custo da 
        %diferença desses histogramas
    end
end

function cost = costh(pointh1, pointh2)
cost = 0;
[fat, ndiv] = size(pointh1); %o número de fatias em que foi dividido o círculo
%centrado em cada ponto (fat) e o número de divisões de cada fatia (ndiv).
for k = 1:fat
    for l = 1:ndiv
        cost_aux = pointh1(k, l) + pointh2(k, l);
        if cost_aux ~= 0
            cost = cost + ((pointh1(k, l) - pointh2(k, l))^2 / cost_aux);
        end
    end
end
cost = cost/2;