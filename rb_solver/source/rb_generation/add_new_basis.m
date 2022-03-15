function [rb] = add_new_basis(rb, mu_idx)

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
X_N = rb.X_N;
N   = rb.N;

mu     = rb.mu;
mu_new = rb.P_train{mu_idx};

B     = rb.pc_data.B;
f     = rb.pc_data.f;
g     = rb.pc_data.g;
theta = rb.DAE.theta;
Y     = rb.pc_data.Y;
R     = rb.R;

solver = rb.options.Solver;

log_start(['adding snapshot nr. ' num2str(N+1)]);

%% Calculate detailed solution
log_middle('assamble right hand side...');
b = f*theta(mu_new) + g;

log_middle('calculate detaild solution...');
N = N + 1;
X_N(:,N) = solver_detailed(B, b, solver);
mu(N) = mu_idx;

%% Calculate riesz representation
log_middle('calculate riesz representation...');
R(:,end+1) = solver_detailed(Y,B.eval(X_N(:,N)),solver);

%% Orthonormalization for stability
log_middle('orthonormalizing the basis functions...');
[X_N, R] = GramSchmidt(X_N,R,rb);

%% Generate new reduced problem
log_middle('construct reduced problem...');
B_ = B.eval2(X_N', X_N);
[BQ, BR] = qr(B_);
B_N.B = B_;
B_N.Q = BQ;
B_N.R = BR;
f_N = X_N'*f;
g_N = X_N'*g;
R_N = Y.inner_product(R, R);

%% Add new basis
rb.X_N = X_N;
rb.B_N = B_N;
rb.f_N = f_N;
rb.g_N = g_N;
rb.R_N = R_N;
rb.R   = R;
rb.N   = N;
rb.mu  = mu;

log_end();
end

