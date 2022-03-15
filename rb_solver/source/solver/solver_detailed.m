function [ x ] = solver_detailed(A,b,solver,local_log_mode)

% ----------------------------------------------------------------------
% REFERENCE:
%  E.Beurer, M.Feuerle, N.Reich, K.Urban
%  "An ultraweak variational method for parmeterized linear 
%  differential-algebraic equations"
%  Ulm University, 2022
%  https://doi.org/10.48550/arXiv.2202.12834
%  https://github.com/mfeuerle/Ultraweak_PDAE
% ----------------------------------------------------------------------

global LOG_MODE;

if nargin == 4
   tmp = LOG_MODE;
   LOG_MODE = local_log_mode;
end

log_start('solving linear system...');

switch upper(solver.Name)
    case 'DIRECT'
        x = A.full()\b;
        flag = -1;
        relres = [];
        iter = [];
    case 'PCG'
        [x,flag,relres,iter] = pcg(@(x) A.eval(x),b, solver.Options{:});
    case 'GMRES'
        [x,flag,relres,iter] = gmres(@(x) A.eval(x),b, solver.Options{:});
    otherwise
        error('Unknown solver.');
end

if flag > 0
    log_warning(sprintf(['solver did not converge. (iteration ' num2str(iter) ', rel. error %e)'],relres));
elseif flag > -1
    log_middle(sprintf(['Iteration ' num2str(iter) ', rel. error %e.'],relres));
end

log_end()

if nargin == 4
    LOG_MODE = tmp;
end
end

