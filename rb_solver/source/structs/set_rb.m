function rb = set_rb(varargin)
    
% ----------------------------------------------------------------------
% REFERENCE:
%  E.Beurer, M.Feuerle, N.Reich, K.Urban
%  "An ultraweak variational method for parmeterized linear 
%  differential-algebraic equations"
%  Ulm University, 2022
%  https://doi.org/10.48550/arXiv.2202.12834
%  https://github.com/mfeuerle/Ultraweak_PDAE
% ----------------------------------------------------------------------

Names = [
    'X_N    '
    'B_N    '
    'f_N    '
    'g_N    '
    'R_N    '
    'R      '
    'N      '
    'mu     '
    'err_N  '
    'P_train'
    'DAE    '
    'options'
    'pc_data'
    ];

stack = dbstack;
% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
    help(stack(1).name);
    return;
end

rb = struct_handler(Names, varargin);
description =@() help(stack(1).name);
rb.help =@() description();

end

