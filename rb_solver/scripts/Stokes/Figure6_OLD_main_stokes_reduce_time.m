% ----------------------------------------------------------------------
% REFERENCE:
%  E.Beurer, M.Feuerle, N.Reich, K.Urban
%  "An ultraweak variational method for parmeterized linear 
%  differential-algebraic equations"
%  Ulm University, 2022
%  https://doi.org/10.48550/arXiv.2202.12834
%  https://github.com/mfeuerle/Ultraweak_PDAE
% ----------------------------------------------------------------------

clear;
close all;
clc;

oldDir = cd('../../..');
initialise_rb_solver;
cd(oldDir);

initialise_example('Stokes');

global LOG_MODE
global TIME_RB_START
LOG_MODE = 0;
TIME_RB_START = tic;
start_time = tic;

save_results = 0;

%% Dimensions
K  = 20;
nx = 15;
ny = 15;
m = 1;

%resolutions = 1:5:K+1;
resolutions = 1:K+1;

%% RBM Options
HomogenMode = 'sigma0'; % Way how the homogenization therm is applied

solver  = set_solver('direct'); % solver for the LGS
% solver  = set_solver('pcg');
% solver  = set_solver('gmres');
% solver  = set_solver('pcg', 0.01*tol, 2*(K+1)*n);

options = set_rb_options('K',           K,...
                         'HomogenMode', HomogenMode,...
                         'Solver',      solver);
                     
%% Time Discretisation
t0 = 0; %Startpunkt
T  = t0 + 1; %Endpunkt
I = [t0 T]; %Intervall
t = linspace(I(1),I(2), K+1);

%% Get Stokes problem
global mu_choice
mu_choice = 3;
global pmax
pmax = ceil(m/2);

%% Either use Full Order system

[E, A, B, ~, v1xx, v2xx, pxx] = stokes_ind2(m, 1, nx, ny);

v10 = - (v1xx(:,2)-0.5);
v20 = (v2xx(:,1)-0.5);
p0  = zeros(size(pxx,1),1);

x0 = [v10; v20; p0];
% x0 = zeros(size(A,1),1);

%% Or use BT reduced system

% load('bt_system_nx15_ny15_m1');
% x0 = ROM.TR\x0_r;
% A = ROM.A;
% E = ROM.E;
% B = ROM.B;

%% Build affine DAE
F = kron(speye(K+1),B);
theta =@(u) u;

DAE = set_DAE('E',E, 'A',A, 'F',F, 'theta',theta, 'x0',x0, 'I',I);

%% RBM generation for different time resolutions of the control
P_validation = eye(K+1,K+1);
err = zeros(length(resolutions),K+1);
err_rel = zeros(length(resolutions),K+1);
err_est = zeros(length(resolutions),K+1);
err_est_rel = zeros(length(resolutions),K+1);

time_pc = tic;
pc_data = precalculate_data(DAE, options);
rb_init = initialising_rb(DAE, pc_data, [], options);
time_pc = toc(time_pc);

time_valid = tic;
sol = cell(K+1,1);
norm_sol = zeros(K+1,1);
for j = 1:K+1
    sol{j} = solver_DAE(DAE, P_validation(:,j), options, pc_data);
    norm_sol(j) = pc_data.X.norm(sol{j}.x_vec);
end
time_valid = toc(time_valid);

times_build = zeros(length(resolutions),1,'uint64');
times_estim = zeros(length(resolutions),1,'uint64');
for i = 1:length(resolutions)
    P_train = discrete_piecewise_linear_basis(resolutions(i),t);
    rb_init.P_train = P_train;
    
    times_build(i) = tic;
    rb = rb_generation_simple(DAE, P_train, options, rb_init);
    times_build(i) = toc(times_build(i));
    
    times_estim(i) = tic;
    for j = 1:K+1
        sol_rb = solver_rb(rb, P_validation(:,j));
        err(i,j) = error_rb(sol_rb, sol{j});
        err_est(i,j) = error_estimate_rb(sol_rb);
        err_rel(i,j) = err(i,j) / norm_sol(j);
        err_est_rel(i,j) = err_est(i,j) / norm_sol(j);
    end
    times_estim(i) = toc(times_estim(i));
    
    fprintf('\n Round %4u; runtime: %010.2f\n\n',i,toc(start_time));
end

%% Handle Results
if save_results == 1
    save(['K' num2str(K) '_n' num2str(size(A,1)) '.mat'], 'err', 'err_est', 'err_est_rel', 'err_rel', 'K', 'resolutions', 'nx', 'ny', 'time_pc', 'time_valid', 'times_build', 'times_estim', 'solver', 'I');
end
plot_time_reduction(resolutions,err_rel,err_est_rel);