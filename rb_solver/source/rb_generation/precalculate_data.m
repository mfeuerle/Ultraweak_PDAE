function [pc_data] = precalculate_data(DAE, options)

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
E  = DAE.E;
A  = DAE.A;
F  = DAE.F;
I  = DAE.I;
x0 = DAE.x0;
n  = DAE.n;
m  = DAE.m;

K           = options.K;
HomogenMode = options.HomogenMode;

log_start('pre-calculating data...');

%% Null space
log_middle('calculating kernel of E^T...');
V = sparse(null(E'));
d = size(V,2);
% V = eye(n);
% d = n;

cN = n*K+d;

%% Matrix products
log_middle('calculating matrix products...');
EE = E*E';
EA = E*A';
AA = A*A';
AAV = AA*V;
EAV = EA*V;
VAAV = V'*AAV;

%% Time Matrices
log_middle('assemble time matrices...');
delta_t = (I(2)-I(1))/K;
t = linspace(I(1),I(2), K+1)';
e = ones(K+1,1);

K_dt = [-e, 2*e, -e];
K_dt = spdiags(K_dt, [-1 0 1], K+1, K+1);
K_dt(1,1) = 1;
K_dt(end,end) = 1;
K_dt = 1/delta_t * K_dt;

L_dt = [e, 4*e, e];
L_dt = spdiags(L_dt, [-1 0 1], K+1, K+1);
L_dt(1,1) = 2;
L_dt(end,end) = 2;
L_dt = delta_t/6 * L_dt;

O_dt = [e, -e];
O_dt = spdiags(O_dt, [-1 1], K+1, K+1);
O_dt(1,1) = -1;
O_dt(end,end) = 1;
O_dt = 0.5 * O_dt;

% M_dt = [e, e];
% M_dt = spdiags(M_dt, [-1 0], K+1, K);
% M_dt = 0.5 * delta_t * M_dt';
% 
% N_dt = [-e, e];
% N_dt = spdiags(N_dt, [0 -1], K+1, K);

%% Bilinear form
% B wird transponiert zur dokumentation aufgestellt
log_middle('assemble bilinear form...');
MT = cell(2,2);
MS = cell(2,2);

MT{1,1} = cell(4,1);
MS{1,1} = cell(4,1);
MT{1,1}{1} = K_dt(1:end-1,1:end-1)';
MT{1,1}{2} = O_dt(1:end-1,1:end-1)';
MT{1,1}{3} = O_dt(1:end-1,1:end-1);
MT{1,1}{4} = L_dt(1:end-1,1:end-1)';
MS{1,1}{1} = EE;
MS{1,1}{2} = EA';
MS{1,1}{3} = EA;
MS{1,1}{4} = AA;

MT{1,2} = cell(2,1);
MS{1,2} = cell(2,1);
MT{1,2}{1} = O_dt(1:end-1,end);
MT{1,2}{2} = L_dt(end,1:end-1)';
MS{1,2}{1} = EAV;
MS{1,2}{2} = AAV;

MT{2,1} = cell(2,1);
MS{2,1} = cell(2,1);
MT{2,1}{1} = O_dt(1:end-1,end)';
MT{2,1}{2} = L_dt(1:end-1,end)';
MS{2,1}{1} = EAV';
MS{2,1}{2} = AAV';

MT{2,2} = cell(1,1);
MS{2,2} = cell(1,1);
MT{2,2}{1} = L_dt(end,end)';
MS{2,2}{1} = VAAV;

B = assemble_matrix(MS, MT);

%% Inner products
% Inner products are identical with the bilinear form
X = B;
Y = B;

%% Linear forms
log_middle('assemble linear forms...');
if isa(F,'function_handle')
    F = zeros(n*(K+1),m);
    for k = 1:K+1
        F(1+n*(k-1):n*k,:) = F(t(k));
    end
end
    
f1 = kron(L_dt(:,1:end-1)',speye(n))*F;
f2 = kron(L_dt(:,end)',V')*F;
f = [f1;f2];

switch HomogenMode
    case 'sigma0'
        Ex0 = E*x0;
        Ax0 = A*x0;
        g1 = kron(L_dt(1,1:end-1)',Ax0) - kron(O_dt(1,1:end-1)',Ex0);
        g2 = L_dt(1,end)*V'*Ax0; % - O_dt(1,end)*V'*Ex0 ;
    case 'const'
        Ax0 = A*x0;
        g1 = kron(delta_t*ones(K,1),Ax0);
        g2 = delta_t*V'*Ax0;
    otherwise
        error('Unknown homogenization function.');
end
g = [g1; g2];


%% Put all together
pc_data = set_pc_data(...
    'B',        B,...
    'f',        f,...
    'g',        g,...
    'cN',       cN,...
    'V',        V,...
    'd',        d,...
    'X',        X,...
    'Y',        Y,...
    'delta_t',  delta_t,...
    't',        t,...
    'K_dt',     K_dt,...
    'O_dt',     O_dt,...
    'L_dt',     L_dt);

log_end();
end




