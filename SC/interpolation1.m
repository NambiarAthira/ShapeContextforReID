

function [b_fR,perimeter]  = interpolation1(b_f,scale,M,str)

% contornos abertos

% if (strcmp(str,'scale'))
%     b_fR     = open_resample(b_f, scale, M, 'scale');
% else
%     b_fR     = open_resample(b_f, scale, M, '');
% end

% contornos fechados

if (strcmp(str,'scale'))
    [b_fR,perimeter]   = open_resample(b_f, scale, M, 'scale');
else
    [b_fR,perimeter]   = open_resample(b_f, scale, M, '');
end


