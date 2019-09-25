%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desenha o "shape context" sobre o ponto em que o utilizador clica
% imag é a imagem inicial
% contourp é o vector de pontos de contorno
% np é o número de pontos do contorno
% radius é o raio do "shape context"
% nslice é o número de fatias em que é dividido o "shape context"
% ndiv é o número de divisões de cada fatia
% type é a forma como cada fatia é dividida (uniforme ou logaritmicamente)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_shapeContext(contourp,np,radius,nslice,ndiv,type)
%hold on,
button = 1;
buttonPressed = 0;

while (button ==1)
    [x,y,button]=ginput(1);

    if button == 1 %add point
        if(buttonPressed == 1)
            undo_plot(gcf, ndiv + nslice);
        end            
        ind = find_closestpoint(contourp, np, x, y);
        hold on,
        draw_circles(contourp(ind), radius, ndiv, type);
        draw_lines(contourp(ind), radius, nslice);
        drawnow
        hold off;
        buttonPressed = 1;
    end
end

%Desenha vários círculos centrados em "point", com raio máximo "radius" e
%dividido "ndiv" vezes uniforme ou logaritmicamente (dependendo do "type")
function draw_circles(point, radius, ndiv, type)
x = real(point);
y = imag(point);
if type == 'u'
    ratio = radius/ndiv;
    for k = 1:ndiv
        draw_circle(x, y, k*ratio);
    end
else
    divider = (2^ndiv)-1;
    for k = 1:ndiv
        draw_circle(x, y, ((2^k-1)/divider)*radius);
    end
end

%Desenha as linhas que vão dividir os círculos, centrados em "point",
%de raio máximo "radius", nas diferentes fatias ("nslice") para formar o 
%"shape context".
function draw_lines(point, radius, nslice)
x1 = real(point);
y1 = imag(point);
ratio = (2*pi)/nslice;
for k = 1:nslice
    x2 = x1+cos(k*ratio)*radius;
    y2 = y1+sin(k*ratio)*radius;
    line([x1 x2],[y1 y2],'Marker','.','LineStyle','-');
end
