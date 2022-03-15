function [ err ] = error_estimate_rb(sol_rb)

% ----------------------------------------------------------------------
% REFERENCE:
%  E.Beurer, M.Feuerle, N.Reich, K.Urban
%  "An ultraweak variational method for parmeterized linear 
%  differential-algebraic equations"
%  Ulm University, 2022
%  https://doi.org/10.48550/arXiv.2202.12834
%  https://github.com/mfeuerle/Ultraweak_PDAE
% ----------------------------------------------------------------------

R_N = sol_rb.rb.R_N;
x_N = sol_rb.x_N;
theta = sol_rb.theta;

r = [theta; 1; -x_N];

err = sqrt(abs(r'*R_N*r));

end

