%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Escreve o índice do ponto à esquerda do ponto
% points é o array com os diferentes pontos
% np é o número de pontos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function label_points(points, np)
for k = 1:np
    s = strcat('-',num2str(k));
    text(real(points(k)),imag(points(k)),s,'HorizontalAlignment','left');
end

