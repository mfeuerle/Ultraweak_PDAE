function [sol_rb] = solver_rb(rb, mu)

% ----------------------------------------------------------------------
% REFERENCE:
%  E.Beurer, M.Feuerle, N.Reich, K.Urban
%  "An ultraweak variational method for parmeterized linear 
%  differential-algebraic equations"
%  Ulm University, 2022
%  https://doi.org/10.48550/arXiv.2202.12834
%  https://github.com/mfeuerle/Ultraweak_PDAE
% ----------------------------------------------------------------------



Q = rb.B_N.Q;
R = rb.B_N.R;
f_N = rb.f_N;
g_N = rb.g_N;

theta = rb.DAE.theta(mu);

b = f_N*theta + g_N;

x_N = R\(Q'*b);

sol_rb.x_N   = x_N;
sol_rb.mu    = mu;
sol_rb.theta = theta;
sol_rb.rb    = rb;

end

