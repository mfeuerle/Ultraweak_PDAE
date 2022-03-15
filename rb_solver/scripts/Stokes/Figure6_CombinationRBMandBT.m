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
global TIME_RB_START
LOG_MODE = 0; % see rb_solver/source/logging/log_start.m
TIME_RB_START = tic; % to have time-stamps in the log

%% Parameters

% Settings for local computation
K = 20;
% Settings as used in the paper
%K = 300;

%DAE
%WARNING: When these parameters are changed a comparison is not possible
%since the BT system was made with these values and is now loaded via a
%mat-file!
nx = 15;
ny = 15;
m = 1;

%% RBM Options
HomogenMode = 'sigma0'; % how the homogenization term is applied
solver  = set_solver('direct'); % solver for the linear system
options = set_rb_options('K',           K,...
                            'HomogenMode', HomogenMode,...
                            'Solver',      solver,...
                            'Nmax',        K+1);
                     
%% Time Discretisation
I = [0 1]; % interval

%% Get Stokes problem and build affine DAE
global mu_choice
mu_choice = 3;    % param for stokes_ind2.m
global pmax
pmax = ceil(m/2); % param for stokes_ind2.m

DAE = cell(2,1);
theta =@(u) u;

%Full order DAE
[E, A, B, ~, v1xx, v2xx, pxx] = stokes_ind2(m, 1, nx, ny);

v10 = - (v1xx(:,2)-0.5);
v20 = (v2xx(:,1)-0.5);
p0  = zeros(size(pxx,1),1);
x0 = [v10; v20; p0];

F = kron(speye(K+1),B);
DAE{1} = set_DAE('E',E, 'A',A, 'F',F, 'theta',theta, 'x0',x0, 'I',I);

%DAE after BT reduction
%BT was done with Matlab-Code by MESS.
%See https://www.mpi-magdeburg.mpg.de/projects/mess for details
load('bt_system_nx15_ny15_m1');
x0 = ROM.TR\x0_r;
A = ROM.A;
E = ROM.E;
B = ROM.B;

F = kron(speye(K+1),B);
DAE{2} = set_DAE('E',E, 'A',A, 'F',F, 'theta',theta, 'x0',x0, 'I',I);

%% Create training set
P_train = rand(K+1,ceil(K*1.5));
P_train = mat2cell(P_train, size(P_train,1), ones(size(P_train,2),1));

%% RB generation through weak greedy algorithm
rb = cell(2,1);
for i=1:2
    rb{i} = rb_generation(DAE{i}, P_train, options);
end

%% Plot
figure
semilogy(rb{1}.err_N,'DisplayName','max rel. err.')
hold on
semilogy(rb{2}.err_N,'DisplayName','max rel. err. with BT')
title('Maximal greedy training error')
xlabel('RB dimension N')
ylabel('max training error')
legend()
