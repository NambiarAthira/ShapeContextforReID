%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Elimina mínimos repetidos de uma matriz de custo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res_cm = elim_min(cm)
res_cm = cm;
dim = size(cm, 1);
for k = 1:dim   %eliminar os mínimos repetidos em cada linha
    minim = min(cm(k,:));
    [row, col] = find(cm(k,:)==minim);
    n_el = numel(row);
    if n_el > 1
        for l = 2:n_el
            res_cm(row(l), col(l)) = res_cm(row(l), col(l)) + 1;
        end
    end
end

for k = 1:dim   %eliminar os mínimos repetidos em cada coluna
    minim = min(cm(:,k));
    [row, col] = find(cm(:,k)==minim);
    n_el = numel(col);
    if n_el > 1
        for l = 2:n_el
            res_cm(row(l), col(l)) = res_cm(row(l), col(l)) + 1;
        end
    end
end