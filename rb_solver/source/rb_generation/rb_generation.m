function [rb] = rb_generation(DAE, P_train, options)

% ----------------------------------------------------------------------
% REFERENCE:
%  E.Beurer, M.Feuerle, N.Reich, K.Urban
%  "An ultraweak variational method for parmeterized linear 
%  differential-algebraic equations"
%  Ulm University, 2022
%  https://doi.org/10.48550/arXiv.2202.12834
%  https://github.com/mfeuerle/Ultraweak_PDAE
% ----------------------------------------------------------------------

global TIME_RB_START;
global LVL_COUNT;
TIME_RB_START  = tic;
LVL_COUNT = [];
warning ('on','all');

Nmax = options.Nmax;
tol  = options.tol;

log_start('constructing a reduced basis...');

%% Pre-calculating and initialising
log_middle('initialising...');

pc_data = precalculate_data(DAE, options);
rb = initialising_rb(DAE, pc_data, P_train, options);

rb = add_new_basis(rb,1);

%% Greedy Iteration
n_train = length(P_train);
mu_left = 2:n_train;
n_left = n_train-1;

while 1
    if n_left <= 0
        log_warning('terminating: out of train space');
        break;
    end
    
    log_middle('calculating error...');
    err = zeros(n_left,1);
    for i = 1:n_left
        sol_rb = solver_rb(rb, P_train{mu_left(i)});
        err(i) = error_estimate_rb(sol_rb); 
    end
    
    [err, pos] = max(err);
    mu_idx = mu_left(pos);
    rb.err_N = [rb.err_N err];
    log_middle(['max. err = ' num2str(err) ' ( at idx = ' num2str(mu_idx) ' )']);
    
    if err < tol
        log_middle(['finished, tol = ' num2str(tol) ' was reached']);
        break;
    end  
    if rb.N >= Nmax
        log_warning('terminating: Nmax reached');
        break;
    end
    
    rb = add_new_basis(rb,mu_idx);
    mu_left(pos) = [];
    n_left = n_left-1;
end
        
log_end('rb generation finished.');
end

