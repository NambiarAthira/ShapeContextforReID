%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Função que elimina n operações de plot anteriores (usada para apagar o
% "shape context" antigo quando se quer analisar um novo.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function undo_plot(imag, n)	
figure(imag);
children = get(gca, 'children');
	
for i = 1 : n
  delete(children(i));
end
