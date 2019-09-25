%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Retorna um array 3D (res) em que cada posi��o desse array � o 
% histograma/matriz com as features (shape contexts) para cada um dos 
% pontos do contorno (contour)
% np is the number of the contour points
% r is the radius od the circle in each contour point
% f is the number of slices in the circle
% n is the number of slices in each "f"
% s is the type of division in "f": 'u' (uniforme) ou 'l' (logar�tmica)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contour=Interp_contour1;
np=M; f=F; n=N2; r=100;
function res = extract_features2(contour,r,np,f,n,s)
res = zeros(f,n,np);
for l = 1:np
    res(:,:,l)=shape_context(l,contour,r,np,f,n,s);
end

 
function sc=shape_context(l,contour,np,r,f,n,s)
sc = zeros(f,n);
x1 = real(contour(l));
y1 = - imag(contour(l)); %necess�rio pois o eixo dos yy est� ao contr�rio do esperado
% v = 1/(np-1); %valor a acrescentar por cada ponto encontrado na respectiva zona
             %de forma a obter um valor normalizado
npoints = 0; %vari�vel para normalizar por n�mero de ocorr�ncias na matriz e 
             %n�o pelo n�mero de pontos
for k = 1:np
    if k~=l  % outro ponto que pertence ao contorno e n�o seja o do �ndice l
        x2 = real(contour(k));
        y2 = - imag(contour(k));
        d = euclid_distance(x1, y1, x2, y2);
        if d <= r   %� sinal que o ponto em k est� dentro do c�rculo
            f_ratio = 2*pi/f;
            ang = atan2(y2-y1, x2-x1);
            if ang < 0
                ang = ang + 2*pi;
            end
            f2add = floor(ang/f_ratio)+1; %qual a fatia a que k pertence 
                             %(1�: 0 a 2*pi/f, exclusiv�; 2�:2*pi/f a
                             %4*pi/f,excl. ...)                                 
            if s == 'u' %divis�o de cada fatia uniformemente
                n_ratio = r/n;
                n2add = floor(d/n_ratio)+1;
            else      %divis�o de cada fatia logaritmicamente
                n_ratio = r/(2^(n)-1);
                n2add = floor(log2((d/n_ratio)+1))+1;                
            end
            n2add = min(n2add, n); %para o caso que d=r, ou est� l� muito 
                                   %perto e o matlab pode aproximar
            f2add = min(f2add, f); %para o caso em que ang seja muito 
                                   %pr�ximo de 0 e depois de somado 2*pi d�
                                   %2*pi devido a aproxima��es do matlab
%            sc(f2add, n2add) = sc(f2add, n2add) + v;
           sc(f2add, n2add) = sc(f2add, n2add) + 1;
            npoints = npoints + 1;
        end
    end
end
end
%if npoints ~= 0
%    sc = sc/npoints;
%end