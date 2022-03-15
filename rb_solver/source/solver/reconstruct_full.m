function [sol] = reconstruct_full(sol_rb)

% ----------------------------------------------------------------------
% REFERENCE:
%  E.Beurer, M.Feuerle, N.Reich, K.Urban
%  "An ultraweak variational method for parmeterized linear 
%  differential-algebraic equations"
%  Ulm University, 2022
%  https://doi.org/10.48550/arXiv.2202.12834
%  https://github.com/mfeuerle/Ultraweak_PDAE
% ----------------------------------------------------------------------

x_vec = reconstruct_reduced(sol_rb);
sol  = recunstruct_detailed(x_vec, sol_rb.rb);

end

