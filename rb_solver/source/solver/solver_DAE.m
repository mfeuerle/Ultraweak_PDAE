function [sol] = solver_DAE(DAE, mu, options, pc_data)

% ----------------------------------------------------------------------
% REFERENCE:
%  E.Beurer, M.Feuerle, N.Reich, K.Urban
%  "An ultraweak variational method for parmeterized linear 
%  differential-algebraic equations"
%  Ulm University, 2022
%  https://doi.org/10.48550/arXiv.2202.12834
%  https://github.com/mfeuerle/Ultraweak_PDAE
% ----------------------------------------------------------------------

%% Gather data
if nargin == 3
    pc_data = precalculate_data(DAE, options);
end
solver = options.Solver;

B       = pc_data.B;
f       = pc_data.f;
g       = pc_data.g;
theta   = DAE.theta;

%% Right hand side
b = f*theta(mu) + g;

%% Solve LGS
[x_vec] = solver_detailed(B,b,solver);
sol = recunstruct_detailed(x_vec, DAE, pc_data, options);
end

