%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Escreve o �ndice do ponto � esquerda do ponto
% points � o array com os diferentes pontos
% np � o n�mero de pontos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function label_points(points, np)
for k = 1:np
    s = strcat('-',num2str(k));
    text(real(points(k)),imag(points(k)),s,'HorizontalAlignment','left');
end

