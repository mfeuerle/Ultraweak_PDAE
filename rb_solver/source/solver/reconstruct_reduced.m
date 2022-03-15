function [x_vec] = reconstruct_reduced(sol_rb)

% ----------------------------------------------------------------------
% REFERENCE:
%  E.Beurer, M.Feuerle, N.Reich, K.Urban
%  "An ultraweak variational method for parmeterized linear 
%  differential-algebraic equations"
%  Ulm University, 2022
%  https://doi.org/10.48550/arXiv.2202.12834
%  https://github.com/mfeuerle/Ultraweak_PDAE
% ----------------------------------------------------------------------

x_N = sol_rb.x_N;
X_N = sol_rb.rb.X_N;

x_vec = X_N * x_N;

end

