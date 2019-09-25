%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coloca a mesma "label" (n�mero) em pontos correspondentes de dois
% vectores. De seguida tra�a segmentos de recta entre o ponto do primeiro 
% vector, mais pr�ximo da posi��o do rato e o seu correspondente.
% contour1 � um dos vectores de pontos de contorno
% contour2 � o outro vector de pontos de contorno
% corresp � a matriz l�gica que d� a correspond�ncia dos pontos 
% np � o n�mero de pontos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function show_correspond(contour1, contour2, corresp, np)
new_cont = zeros(np, 1);
for k = 1:np
    pos = find(corresp(k,:));
    new_cont(k) = contour2(pos);
end

hold on,
label_points(contour1, np);
label_points(new_cont, np);
hold off
draw_corresp(contour1, new_cont, np);

%Desenha segmentos de recta entre pontos correspondentes tendo em conta a
%posi��o do rato.
function draw_corresp(cont1, cont2, np)
button = 1;
ind_aux = 0;

while (button == 1)
    [x,y,button]=ginput(1);
            
    ind = find_closestpoint(cont1, np, x, y);
    if(ind ~= ind_aux)
        if(ind_aux ~= 0 )
            undo_plot(gcf, 1);
        end
        hold on,
        line([real(cont1(ind)) real(cont2(ind))],[imag(cont1(ind)) imag(cont2(ind))],'Marker','.','LineStyle','-','Color','black','LineWidth',1);
        drawnow
        hold off;
        ind_aux = ind;
    end
end

% %FICAVA MUITO CONFUSO
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Desenha linhas a preto a unir pontos correspondentes
% % contour1 � um dos vectores de pontos de contorno
% % contour2 � o outro vector de pontos de contorno
% % corresp � a matriz l�gica que d� a correspond�ncia dos pontos 
% % dim � o n�mero de pontos
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function show_correspond(contour1, contour2, corresp, dim)
% hold on
% for k = 1:dim
%     pos = find(corresp(k,:));
%     line([real(contour1(k)) real(contour2(pos))],[imag(contour1(k)) imag(contour2(pos))],'Marker','.','LineStyle','-','Color','black','LineWidth',1);
% end
% hold off