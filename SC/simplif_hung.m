%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Retorna uma matriz de atribui��es e o seu custo tendo em conta a matriz 
% de custo "cm".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [atrib_mat, cost] = simplif_hung(cm)
dim = size(cm, 1); %vamos assumir que a matriz � quadrada (mesmo n�mero de 
                   %"shape contexts" que geram a matriz de custo
atrib_mat(1:dim, 1:dim) = false;
n_atrib = 0;
cost = 0;
%a verifica��o em rela��o a um threshold � para todas as posi��es da matriz
%ou para os m�nimos?
while n_atrib ~= dim
    c_aux = find_min_lines(cm, dim);
    c_aux = c_aux + find_min_cols(cm, dim);
    [n, c, at_mat, cost_aux] = find_twos(c_aux, cm, atrib_mat, dim);
    cm = c;
    atrib_mat = at_mat;
    cost = cost + cost_aux;
    n_atrib = n_atrib + n;
    if n_atrib == dim - 1
        [line, col, cost_aux] = final_atrib(cm, dim);
        atrib_mat(line, col) = true;
        cost = cost + cost_aux;
        n_atrib = dim;
    else
        if n == 0
            disp('Could not found any new attributions');
            return;
        end
    end    
end

% %This function places the number 1 in the position(s) with the minimum value
% %for each line of the matrix "mat", with dimensions "dimxdim".
% function mlin = find_min_lines(mat, dim)
% mlin = zeros(dim);
% min = inf;
% %pos = 0;
% for k = 1:dim
%     for l = 1:dim
% %        if mat(k, l) ~= inf
%          if mat(k, l) < min
%              min = mat(k, l);
% %            pos = l;
%          end
% %        end
%     end
%     if min ~= inf   %encontrou n�meros sem ser inf na linha k
%         for l = 1:dim  %� necess�rio este ciclo pois na mesma linha pode 
%                        %haver mais que um m�nimo
%             if mat(k, l) == min
%                 mlin(k, l) = 1;
%             end
%         end
%         min = inf;
% %        pos = 0;
%     end
% end

%Version using Matlab capabilities
%This function places the number 1 in the position(s) with the minimum value
%for each line of the matrix "mat", with dimensions "dimxdim".
function mlin = find_min_lines(mat, dim)
mlin(1:dim, 1:dim) = false;
for k = 1:dim
    minim = min(mat(k,:));
    if minim ~= inf   %encontrou n�meros sem ser inf na linha k
        mlin(k,:) = (mat(k,:) == minim); %coloca 1s na linha k nas posi��es 
                                       %com valores igual ao m�nimo
    end
end

% %This function places the number 1 in the position(s) with the minimum value
% %for each column of the matrix "mat", with dimensions "dimxdim".
% function mcol = find_min_cols(mat, dim)
% mcol = zeros(dim);
% min = inf;
% %pos = 0;
% for k = 1:dim
%     for l = 1:dim
% %        if mat(k, l) ~= inf
%         if mat(l, k) < min
%             min = mat(l, k);
% %            pos = l;
%         end
% %        end
%     end
%     if min ~= inf   %encontrou n�meros sem ser inf na coluna k
%         for l = 1:dim  %� necess�rio este ciclo pois na mesma coluna pode 
%                        %haver mais que um m�nimo
%             if mat(l, k) == min
%                 mcol(l, k) = 1;
%             end
%         end
%         min = inf;
% %        pos = 0;
%     end
% end

%Version using Matlab capabilities
%This function places the number 1 in the position(s) with the minimum value
%for each column of the matrix "mat", with dimensions "dimxdim".
function mcol = find_min_cols(mat, dim)
mcol(1:dim, 1:dim) = false;
for k = 1:dim
    minim = min(mat(:,k));
    if minim ~= inf   %encontrou n�meros sem ser inf na coluna k
        mcol(:,k) = (mat(:,k) == minim); %coloca 1s na coluna k nas posi��es 
                                       %com valores igual ao m�nimo
    end
end

%For each "2" this function finds in "two_mat" it replaces the elements of the 
%corresponding line and column in "cost_mat" with inf and places a "1" in the 
%at_mat. It also returns the cost of the elements selected, "cost", and 
%its number "n"   
function [n, c_mat, at_mat, cost] = find_twos(two_mat, cost_mat, atr_mat, dim)
n = 0;
cost = 0;
at_mat = atr_mat;
c_mat = cost_mat;
for k = 1:dim
    l = 1;
    while l ~= dim + 1
        if two_mat(k, l) == 2 && cost_mat(k, l) ~= inf %n�o foi eliminado numa itera��o anterior
            cost = cost + cost_mat(k, l);
            n = n + 1;
            at_mat(k, l) = true;
            c_mat = place_inf(cost_mat, k, l, dim);
            cost_mat = c_mat; %we update the matrix and keep searching it for 2s
            l = dim;
        end
        l = l +1;
    end
end

%This function places "inf" in line "line" and column "col".
function res_mat = place_inf(mat, line, col, dim)
for k = 1:dim
    mat(line, k) = inf;
    mat(k, col) = inf;
end
res_mat = mat;

%This function returns the "line", column ("col") and "cost" of the single 
%element in "cost_m" different than inf.
function [line, col, cost] = final_atrib(cost_m, dim)
for k = 1:dim
    for l = 1:dim
        if cost_m(k, l) ~= inf
            line = k;
            col = l;
            cost = cost_m(k, l);
        end
    end
end