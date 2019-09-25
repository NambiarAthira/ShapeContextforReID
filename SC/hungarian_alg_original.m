%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Retorna uma matriz de atribuições e o seu custo tendo em conta a matriz 
% de custo "cm". (baseado nos pdfs-hungarian.pdf e munkres-calcular minimo... 
%linhas para cobrir zeros.pdf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Returns an array of tasks and their cost taking into account the matrix 
% Cost "cm". (based on hungarian.pdf pdfs-and minimum-munkres calculate ... 
% lines to cover zeros.pdf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [atrib_mat, cost] = hungarian_alg(cm)
dim = size(cm, 1); %vamos assumir que a matriz é quadrada (mesmo número de 
                   %"shape contexts" que geram a matriz de custo
cm_aux = cm;

%STEP 1-subtrair o mínimo de todas as linhas
for k = 1:dim
    minim = min(cm_aux(k,:));
    cm_aux(k,:) = cm_aux(k,:) - minim;
end

%STEP 2-subtrair o mínimo de todas as colunas
for k = 1:dim
    minim = min(cm_aux(:,k));
    cm_aux(:,k) = cm_aux(:,k) - minim;
end

%STEP 3-Cobrir todos os zeros de cm com o mínimo número de linhas e colunas
%até esse mínimo ser igual ao número de linhas (e colunas) da matriz cm
[cov, stZ] = cover_zeros(cm_aux, dim);
%k = n_cov(cov, dim); %este número talvez seja o número de elementos em
%starZ. Segundo o artigo é.
k = numel(find(stZ == true));
prZ(1:dim, 1:dim) = false;
while (k ~= dim)
    cov_aux = cov;
    stZ_aux = stZ;
    prZ_aux = prZ;
    %STEP 3.1
    [cov, stZ, prZ] = step1_cov(stZ_aux, prZ_aux, cov_aux, cm_aux, dim);
    %k = n_cov(cov, dim);
    k = numel(find(stZ == true));
    if k ~= dim
        min_unc_numb = min_unc(cm_aux, cov, dim);
        cm_aux2 = change_mat(cm_aux, cov, min_unc_numb, dim);
        cm_aux = cm_aux2;
    end
end

%STEP 4-Fazer as atribuições
[atrib_mat, cost] = atrib(cm_aux, cm, dim);

function [cov, starZ] = cover_zeros(mat, dim)
cov = zeros(dim);
starZ(1:dim, 1:dim) = false;

[row, col] = find(mat == 0);
s = size(row, 1);
for k = 1:s
    if numel(find(starZ(row(k),:))) == 0 && numel(find(starZ(:,col(k)))) == 0
        starZ(row(k), col(k)) = true;
    end
end

[row, col] = find(starZ);
if numel(col) ~= 0
    cov(:,col) = cov(:,col) + 1;
end

function [res_cov, res_starZ, res_primeZ] = step1_cov(starZ, primeZ, cov, mat, dim)
prim_Z = primeZ;
while(true)
    c_aux = mat + cov;
    [row, col] = find(c_aux == 0);
    if numel(row) == 0
        res_cov = cov;
        res_starZ = starZ;
        res_primeZ = prim_Z;
        return;
    else
        prim_Z(row(1), col(1)) = true;
        [row2, col2] = find(starZ(row(1),:));
        if numel(col2) == 0
            %STEP 3.2
            [res_cov, res_starZ, res_primeZ] = step2_cov(starZ, prim_Z, row(1), col(1), cov, mat, dim);
            cov = res_cov; %é necessário actualizar a matriz de cobertura
            starZ = res_starZ;
            prim_Z = res_primeZ;
        else
            cov(row(1),:) = cov(row(1),:) + 1;
            cov(:,col2) = cov(:,col2) - 1;
        end
    end
end

function [res_cov, res_starZ, res_primeZ] = step2_cov(starZ, primeZ, r, c, cov, mat, dim)
aux_starZ = starZ;
seq_star = [];
n_seq = 1;
seq_prim(n_seq) = r + 1i*c;

[rs, cs] = find(starZ(:, c));
cs = c; %o find não devolve a coluna certa
while numel(rs) ~= 0
    seq_star(n_seq) = rs + 1i*cs;
    [rp, cp] = find(primeZ(rs,:));
    rp = rs; %o find não devolve a linha certa
%     if numel(cp) == 0  %segundo o texto em que te baseaste esta verificação deve ser inútil
%         disp('Error1: hungarian_alg');
%         return;
%     else
    n_seq = n_seq + 1;
    seq_prim(n_seq) = rp + 1i*cp;
%     end
    [rs, cs] = find(starZ(:, cp));
    cs = cp;
end

%Unstar each zero of the star sequence
n_seq = numel(seq_star);
if n_seq ~= 0
    for k = 1:n_seq
        aux_starZ(real(seq_star(k)), imag(seq_star(k))) = false;
    end
end

%Star each zero of the prime sequence
n_seq = numel(seq_prim);
for k = 1:n_seq
    aux_starZ(real(seq_prim(k)), imag(seq_prim(k))) = true;
end

%Erase all primes
aux_primeZ(1:dim, 1:dim) = false;

cov_chang = cov;
%Uncover every row
for k = 1:dim
    row_cov = find(cov_chang(k,:));
    if numel(row_cov) == dim
        cov_chang(k,:) = cov_chang(k,:) - 1;
    end
end

%Cover every column containing a star zero
[rz, cz] = find(aux_starZ); %há pelo menos um elemento
nz = numel(cz);
for k = 1:nz  %vou ter de confiar que não há mais que um elemento a 1 em starZ na mesma coluna
    if numel(find(cov_chang(:,cz(k)))) ~= dim %a coluna ainda não está coberta
        cov_chang(:,cz(k)) = cov_chang(:,cz(k)) + 1;
    end
end

%STEP 3.1
[res_cov, res_starZ, res_primeZ] = step1_cov(aux_starZ, aux_primeZ, cov_chang, mat, dim);

%Calcula o número de linhas e colunas cobertas
function num_cov = n_cov(cov, dim)
num_cov = 0;
cov_aux = cov;
for k = 1:dim
    if numel(find(cov_aux(k,:))) == dim
        num_cov = num_cov + 1;
        cov_aux(k,:) = inf;
    end
end

for k = 1:dim
    if numel(find(cov_aux(:,k))) == dim
        if numel(find(cov_aux(:,k) == inf)) ~= dim %a coluna não foi já coberta por linhas que já foram contadas.
            num_cov = num_cov + 1;
        else %a matriz está toda coberta
            return;
        end
    end
end

%Descobre o elemento mínimo de cm que não está coberto por cov
function minim_uncov = min_unc(cm, cov, dim)
minim_uncov = inf;
for k = 1:dim
    for l = 1:dim
        if cov(k,l) == 0
            aux = cm(k,l);
            if aux < minim_uncov
                minim_uncov = aux;
            end
        end
    end
end

%Altera cm subtraindo minim_uncov de cada número que não está coberto e
%somando a mesma quantia a cada número que se encontra coberto por linha e
%coluna. A cobertura é avaliada de acordo com co
function res_cm = change_mat(cm, cov, minim_uncov, dim)
res_cm = cm;
for k = 1:dim
    for l = 1:dim
        if cov(k,l) == 0 %não está coberto
            res_cm(k,l) = res_cm(k,l) - minim_uncov;
        else
            if cov(k,l) == 2 %coberto por linha e coluna
                res_cm(k,l) = res_cm(k,l) + minim_uncov;
            end
        end
    end
end

%Função que cria a matriz lógica com as correspondências e calcula o custo
%total.
function [atrib_mat, cost] = atrib(cm_chang, cm, dim)
atrib_mat(1:dim, 1:dim) = false;
cost = 0;
n_atrib = 0;

while n_atrib ~= dim
    n_aux = -1;
    while n_aux ~= n_atrib
        n_aux = n_atrib;
        k = 1;
        while (k ~= dim + 1) && (n_atrib ~= dim) %atribuições únicas linha a linha
            [row, col] = find(cm_chang(k,:) == 0);
            if numel(col) == 1
                n_atrib = n_atrib + 1;
                atrib_mat(k,col) = true;
                cost = cost + cm(k,col);
                cm_aux = place_inf(cm_chang, k, col);
                cm_chang = cm_aux;
            end
            k = k + 1;
        end
        if n_atrib ~= dim
            k = 1;
            while (k ~= dim + 1) && (n_atrib ~= dim)%atribuições únicas coluna a coluna
                [row, col] = find(cm_chang(:,k) == 0);
                if numel(row) == 1
                    n_atrib = n_atrib + 1;
                    atrib_mat(row,k) = true;
                    cost = cost + cm(row,k);
                    cm_aux = place_inf(cm_chang, row, k);
                    cm_chang = cm_aux;
                end
                k = k + 1;
            end
        else
            return;
        end
    end
    k = 1;
    while k ~= dim + 1 %atribuições aleatória entre as escolhas possíveis (a primeira que apareça)
        [row, col] = find(cm_chang(k,:) == 0);
        if numel(col) > 0
            n_atrib = n_atrib + 1;
            atrib_mat(k,col(1)) = true;
            cost = cost + cm(k,col(1));
            cm_aux = place_inf(cm_chang, k, col(1));
            cm_chang = cm_aux;
            k = dim;
        end
        k = k + 1;
    end
end

%This function places "inf" in line "line" and column "col".
function res_mat = place_inf(mat, line, col)
res_mat = mat;

res_mat(line, :) = inf;
res_mat(:, col) = inf;