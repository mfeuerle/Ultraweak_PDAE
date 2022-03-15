function [sol] = recunstruct_detailed(x_vec, DAE, pc_data, options)

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
if nargin == 2
    rb = DAE;
    DAE = rb.DAE;
    pc_data = rb.pc_data;
    options = rb.options;
end

E  = DAE.E;
A  = DAE.A;
x0 = DAE.x0;
n  = DAE.n;

K = options.K;
HomogenMode = options.HomogenMode;

V = pc_data.V;
d = pc_data.d;
delta_t = pc_data.delta_t;
t = pc_data.t;

%% Reconstruct
x = [reshape(x_vec(1:n*K),[n,K]),...
    V*reshape(x_vec(n*K+1:end),[d,1])];
x = eval(x);

switch HomogenMode
    case 'sigma0'
        x(:,1,2) = x(:,1,2) + x0;
    case 'const'
        x = x + x0;
    otherwise
        error('Unknown homogenization function.');
end

sol.x = x;
sol.t = t;
sol.x_vec = x_vec;
sol.DAE = DAE;
sol.options = options;
sol.pc_data = pc_data;

    function x_val = eval(x)
        x_val = zeros(size(A,1),K+1,2);
        
        x_val(:,1,1) = nan;
        x_val(:,1,2) = E'*x(:,1)/delta_t - A'*x(:,1) ...
            - E'*x(:,2)/delta_t;
        
        for k = 2:K
            x_val(:,k,1) =   E'*x(:,k-1)/delta_t ...
                - E'*x(:,k)/delta_t - A'*x(:,k);
            x_val(:,k,2) =   E'*x(:,k)/delta_t - A'*x(:,k) ...
                - E'*x(:,k+1)/delta_t;
        end
        
        x_val(:,K+1,1) =   E'*x(:,K)/delta_t ...
            - E'*x(:,K+1)/delta_t - A'*x(:,K+1);
        x_val(:,K+1,2) =   nan;
    end
end


