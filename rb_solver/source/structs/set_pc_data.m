function options = set_pc_data(varargin)
        
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
    'B      '
    'f      '
    'g      '
    'cN     '
    'V      '
    'd      '
    'X      '
    'Y      '
    'delta_t'
    't      '
    'K_dt   '
    'O_dt   '
    'L_dt   '
    ];

stack = dbstack;
% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
    help(stack(1).name);
    return;
end

options = struct_handler(Names, varargin);
description =@() help(stack(1).name);
options.help =@() description();

end

