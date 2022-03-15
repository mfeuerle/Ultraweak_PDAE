function rb = initialising_rb(DAE, pc_data, P_train, options)

% ----------------------------------------------------------------------
% REFERENCE:
%  E.Beurer, M.Feuerle, N.Reich, K.Urban
%  "An ultraweak variational method for parmeterized linear 
%  differential-algebraic equations"
%  Ulm University, 2022
%  https://doi.org/10.48550/arXiv.2202.12834
%  https://github.com/mfeuerle/Ultraweak_PDAE
% ----------------------------------------------------------------------

cN = pc_data.cN;
Y  = pc_data.Y;
f  = pc_data.f;
g  = pc_data.g;

solver = options.Solver;

m  = DAE.m;

log_start('prepair RBM...');

%% Riesz representations for calculating the dual norm of the residuum
R = zeros(cN,m+1);
for i = 1:m     % lässt sich parallelisieren (probleme mit logging)
    log_middle(['calculating riesz representations (' num2str(i) ' of ' num2str(m+1) ') ...']);
    R(:,i) = solver_detailed(Y,f(:,i),solver,2);
end
log_middle(['calculating riesz representations (' num2str(m+1) ' of ' num2str(m+1) ') ...']);
R(:,m+1) = solver_detailed(Y,g,solver,2);


%% Setup rest
N = 0;

B_N.B = zeros(N,N);
B_N.Q = zeros(N,N);
B_N.R = zeros(N,N);

X_N = zeros(cN,N);
f_N = zeros(N,m);
g_N = zeros(N,1);

R_N = Y.inner_product(R,R);

mu = zeros(N,1);

err_N = zeros(1,N);

rb = set_rb('X_N',      X_N,...
            'B_N',      B_N,...
            'f_N',      f_N,...
            'g_N',      g_N,...
            'R_N',      R_N,...
            'R',        R,...
            'N',        N,...
            'mu',       mu,...
            'err_N',    err_N,...
            'P_train',  P_train,...
            'DAE',      DAE,...
            'options',  options,...
            'pc_data',  pc_data);

log_end();
end

