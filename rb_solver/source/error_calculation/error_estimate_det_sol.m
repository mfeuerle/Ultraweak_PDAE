function [err] = error_estimate_det_sol(sol,DAE,B,mu,options,K2)

% ----------------------------------------------------------------------
% REFERENCE:
%  E.Beurer, M.Feuerle, N.Reich, K.Urban
%  "An ultraweak variational method for parmeterized linear 
%  differential-algebraic equations"
%  Ulm University, 2022
%  https://doi.org/10.48550/arXiv.2202.12834
%  https://github.com/mfeuerle/Ultraweak_PDAE
% ----------------------------------------------------------------------

%% Prep for finer time discretization

DAE2 = DAE;
options.K = K2;
t = linspace(DAE.I(1),DAE.I(2), K2+1);
DAE2.F = kron(speye(K2+1),B);
DAE2.theta =@(u) u(t');
pc_data = precalculate_data(DAE2, options);
rhs = pc_data.f*DAE2.theta(mu) + pc_data.g;

%% Space matrices

E = DAE.E;
A = DAE.A;
N = pc_data.V;
V = pc_data.V;

EE = E*E';
AE = A*E';
EA = E*A';
AA = A*A';
NEE = N'*EE;
NAE = N'*AE;
NEA = N'*EA;
NAA = N'*AA;

AAV = AA*V;
EAV = EA*V;
VAAV = V'*AAV;

%% time matrices

K1 = length(sol.t)-1;

t1 = sol.t;
t2 = linspace(DAE.I(1),DAE.I(2),K2+1)';

delta_t1 = t1(2)-t1(1);
delta_t2 = t2(2)-t2(1);

K_dt = sparse(K2+1,K1+1);   %(\dot\sigma_2, \dot\sigma_1)
L_dt = sparse(K2+1,K1+1);   %(\sigma_2, \sigma_1)
O1_dt = sparse(K2+1,K1+1);  %(\sigma_2, \dot\sigma_1)
O2_dt = sparse(K2+1,K1+1);  %(\dot\sigma_2, \sigma_1)

%% Calculation of (struct) B for final linear system of equations

k2_prev = 1;
for k1 = 1:K1
    k2 = k2_prev;
    while t2(k2+1) <= t1(k1)
        k2 = k2+1;
    end
    k2_prev = k2;
    
    while t2(k2) < t1(k1+1)
        a = max(t1(k1),t2(k2));
        b = min(t1(k1+1),t2(k2+1));
        
        x1 = t1(k1)/delta_t1;
        x2 = t2(k2)/delta_t2;
        
        int1 = (b^3-a^3)/(3*delta_t1*delta_t2);
        int2 = (b^2-a^2)/(2*delta_t1);
        int3 = (b^2-a^2)/(2*delta_t2);
        
        L = zeros(2);
        L(1,1) =  (1+x1)*(1+x2)*(b-a) + int1 - (1+x2)*int2 - (1+x1)*int3;
        L(1,2) =     -x1*(1+x2)*(b-a) - int1 + (1+x2)*int2 +     x1*int3;
        L(2,1) = -(1+x1)*   x2 *(b-a) - int1 +     x2*int2 + (1+x1)*int3;
        L(2,2) =      x1*   x2 *(b-a) + int1 -     x2*int2 -     x1*int3;
        L_dt([k2 k2+1],[k1 k1+1]) = L_dt([k2 k2+1],[k1 k1+1]) + L;
        
        K = (b-a)/(delta_t1*delta_t2) * [1, -1; -1, 1];
        K_dt([k2 k2+1],[k1 k1+1]) = K_dt([k2 k2+1],[k1 k1+1]) + K;
        
        O1 = int3/delta_t1 * [1, -1; -1, 1] + (b-a)/delta_t1 * [-1-x2, 1+x2; x2, -x2];
        O1_dt([k2 k2+1],[k1 k1+1]) = O1_dt([k2 k2+1],[k1 k1+1]) + O1;
        
        O2 = int2/delta_t2 * [1, -1; -1, 1] + (b-a)/delta_t2 * [-1-x1, x1; 1+x1, -x1];
        O2_dt([k2 k2+1],[k1 k1+1]) = O2_dt([k2 k2+1],[k1 k1+1]) + O2;
        
        k2 = k2+1;
        if k2 > K2
            break;
        end
    end
end

MT = cell(2,2);
MS = cell(2,2);

MT{1,1} = cell(4,1);
MS{1,1} = cell(4,1);
MT{1,1}{1} = K_dt(1:end-1,1:end-1);
MT{1,1}{2} = O1_dt(1:end-1,1:end-1);
MT{1,1}{3} = O2_dt(1:end-1,1:end-1);
MT{1,1}{4} = L_dt(1:end-1,1:end-1);
MS{1,1}{1} = EE;
MS{1,1}{2} = EA';
MS{1,1}{3} = EA;
MS{1,1}{4} = AA;

MT{1,2} = cell(2,1);
MS{1,2} = cell(2,1);
MT{1,2}{1} = O2_dt(1:end-1,end);%O2_dt(1:end-1,end)';
MT{1,2}{2} = L_dt(1:end-1,end);
MS{1,2}{1} = EAV;
MS{1,2}{2} = AAV;

MT{2,1} = cell(2,1);
MS{2,1} = cell(2,1);
MT{2,1}{1} = O1_dt(end,1:end-1);
MT{2,1}{2} = L_dt(end,1:end-1);
MS{2,1}{1} = EAV';
MS{2,1}{2} = AAV';

MT{2,2} = cell(1,1);
MS{2,2} = cell(1,1);
MT{2,2}{1} = L_dt(end,end);
MS{2,2}{1} = VAAV;

B = assemble_matrix(MS, MT);

%% Calculation of error estimate

res = B.eval(sol.x_vec) - rhs;
err = sqrt(res'*solver_detailed(pc_data.Y,res,options.Solver));

end

