function [L2product] = L2_product_discrete(f,g,t)
% calculates the L2 product
%   (f,g)_L2
% of two diskrete functions f,g, given at the times t with the
% trapezoidal rule.

% ----------------------------------------------------------------------
% REFERENCE:
%  E.Beurer, M.Feuerle, N.Reich, K.Urban
%  "An ultraweak variational method for parmeterized linear 
%  differential-algebraic equations"
%  Ulm University, 2022
%  https://doi.org/10.48550/arXiv.2202.12834
%  https://github.com/mfeuerle/Ultraweak_PDAE
% ----------------------------------------------------------------------

if isa(f,'function_handle')
    f = f(t);
end
if isa(g,'function_handle')
    g = g(t);
end

fg = sum(f.*g,1);
L2product = trapz(t,fg);

end
