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

global LOG_MODE
LOG_MODE = 0;   % see rb_solver/source/logging/log_start.m
global TIME_RB_START
TIME_RB_START = tic; % to have time-stamps in the log

%% RLC Circuit Problem
L = 1;    % inductivity in Henry
C = 1;    % capacity in Farad
R = 1e-8; % resistance in Ohm

E = [1   0   0   0;
     0   1   0   0;
     0   0   0   0;
     0   0   0   0];
A = [0   0 1/L   0;
     1/C 0   0   0;
     R   0   0  -1;
     0   1   1   1];
B = [0;  0;  0; -1];

%% Voltage source and corresponding analytical solution
U_hat = 5;           % voltage amplitude in Volt
omega = 1/sqrt(L*C); % frequency in 1/s
phi_U = 0;           % phase shift at t0 in rad

[x_analy, f_Vs] = get_analytic_solution(R,L,C,U_hat,phi_U,omega);

%% Time intervall and consistent (parameter independent) initial value
t0 = 0; % starting point
T  = t0 + 2*2*pi/omega; % end
I = [t0 T]; % interval

x0 = x_analy(I(1));

%% Number of time steps K
% Settings that should be executable locally in reasonable time:
K  = [10,20,50,100,200,500,1000];
% Settings we used for the paper (on a cluster):
% K = [10,20,50,100,200,500,1000,2000,5000,10000,20000];

%% Run Problem
err = zeros(length(K),1);
err_est = zeros(length(K),1);

for i = 1:length(K)
    %% Assamble affine decomposition of the rhs
    t = linspace(I(1),I(2), K(i)+1);
    F = kron(speye(K(i)+1),B);
    theta =@(u) u(t');

    %% Assamble the full problem
    DAE = set_DAE('E',E, 'A',A, 'F',F, 'theta',theta, 'x0',x0, 'I',I);
    
    %% Calculate Petrov-Galerking approximation
    HomogenMode = 'const'; % Way how the homogenization of x0 is applied
    solver  = set_solver('direct'); % solver for the LGS
    options = set_rb_options('K',K(i), 'HomogenMode',HomogenMode, 'Solver',solver);

    sol = solver_DAE(DAE, f_Vs, options);
    
    %% Calculate error and error estimator
    % for the error estimator, a finer approximation (resolution) of the 
    % infinit spaces X,Y is needed
    res_reference = 20*K(i); 

    err(i) = error_det_sol(sol,x_analy,100);
    err_est(i) = error_estimate_det_sol(sol,DAE,B,f_Vs,options,res_reference);
end

%% Calculate relative L2 error
ZERO.t = linspace(I(1),I(2),max(K)+1);
ZERO.x = zeros(size(E,2),length(ZERO.t));
norm_exact = get_L2_error(ZERO,x_analy,100); % calcs ||0-x_analy||_L2

err = err./norm_exact;
err_est = err_est./norm_exact;

%% Plots
figure()

loglog(K,err,'ks-', 'DisplayName','relative error')
hold on
loglog(K,err_est,'rd--', 'DisplayName','relative error estimator')

title('h-convergence and error estimation')
xlabel('number of time steps K')
ylabel('relative error/estimator')
legend()

