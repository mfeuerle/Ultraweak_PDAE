function [ err ] = error_rb(sol_rb, sol_det)

% ----------------------------------------------------------------------
% REFERENCE:
%  E.Beurer, M.Feuerle, N.Reich, K.Urban
%  "An ultraweak variational method for parmeterized linear 
%  differential-algebraic equations"
%  Ulm University, 2022
%  https://doi.org/10.48550/arXiv.2202.12834
%  https://github.com/mfeuerle/Ultraweak_PDAE
% ----------------------------------------------------------------------

rb      = sol_rb.rb;
theta   = sol_rb.theta;

solver = rb.options.Solver;

B       = rb.pc_data.B;
f       = rb.pc_data.f;
g       = rb.pc_data.g;


%% Calculate detailed solution
if nargin < 2
    b = f*theta + g;
    [x_vec] = solver_detailed(B,b,solver);
else
    x_vec = sol_det.x_vec;
end

%% Fehlerberechnung (Methode 1)
x_vec_rb = reconstruct_reduced(sol_rb);
err = rb.pc_data.X.norm(x_vec - x_vec_rb);

%% Fehlerberechnung (Methode 2)
% K = rb.options.K;
% delta_t = rb.pc_data.delta_t;
% 
% sol_rb = recunstruct_detailed(x_vec_rb, rb);
% sol    = recunstruct_detailed(x_vec, rb);
% 
% diff = sol.x - sol_rb.x;
% 
% err2 = 0;
% for k = 1:K
%     err2 = err2 + ...
%         delta_t * (diff(:,k,2)'*diff(:,k,2) + diff(:,k+1,1)'*diff(:,k+1,1))/2;
% end
% err2 = sqrt(abs(err2));

end

