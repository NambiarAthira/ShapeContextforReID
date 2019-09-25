function Vo = my_resample(Vi,scale,Mo,str)

%function [Vo] = my_resample(Vi)
%Mo = length(Vi);
Mi = length(Vi);
for i=1:Mi-1
    d(i) = sqrt((Vi(i+1,:) - Vi(i,:)) * (Vi(i+1,:) - Vi(i,:))');
end;
d(Mi) = sqrt((Vi(Mi,:) - Vi(1,:)) * (Vi(Mi,:) - Vi(1,:))');
indices=[1:length(d)];
perimeter = sum(d);

if (strcmp(str,'scale'))
    Mo = round(perimeter / scale);
end

if Mo < 6
    disp('Contour has collapsed');
    return;
end;

E_med = perimeter / Mo;

Vo(1,:) = Vi(1,:);

i = 1; Ei = 0;

for n=2:Mo
    D = E_med + Ei;
    j = i;
    aux = 0;
    while aux < D
        aux = aux + d(j);
        j = j + 1;
    end;
    i = j - 1;
    if j == (Mi + 1)
        j = 1;         % curva fechada
        %j = Mi - 1;   % curva aberta 
    end;
    Ei = D - aux + d(i);
    Vo(n,:) = Vi(i,:) + (Vi(j,:) - Vi(i,:)) * Ei / d(i);
end;
