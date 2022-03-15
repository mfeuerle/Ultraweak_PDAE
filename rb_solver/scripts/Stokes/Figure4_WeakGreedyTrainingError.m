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
LOG_MODE = 0;   % see rb_solver/source/logging/log_start.m
TIME_RB_START = tic; % to have time-stamps in the log

%% Parameters

% Settings for local computation
K_min = 5;
% Settings as used in the paper
%K_min = 75;

% Ks should have hierarchical order with factors 1,2,4,...
K = [K_min, 2*K_min, 4*K_min];

%DAE
nx = 15;
ny = 15;
m = 1;

%% RBM Options
HomogenMode = 'sigma0'; % how the homogenization term is applied

solver  = set_solver('direct'); % solver for the linear system

options = cell(length(K),1);
for i=1:length(K)
    options{i} = set_rb_options('K',           K(i),...
                                'HomogenMode', HomogenMode,...
                                'Solver',      solver,...
                                'Nmax',        K(i)+1);
end
                     
%% Time Discretisation
I = [0 1]; % interval

%% Get Stokes problem
global mu_choice
mu_choice = 3;    % param for stokes_ind2.m
global pmax
pmax = ceil(m/2); % param for stokes_ind2.m

[E, A, B, ~, v1xx, v2xx, pxx] = stokes_ind2(m, 1, nx, ny);

v10 = - (v1xx(:,2)-0.5);
v20 = (v2xx(:,1)-0.5);
p0  = zeros(size(pxx,1),1);

x0 = [v10; v20; p0];

%% Build affine DAE
DAE = cell(length(K),1);
for i=1:length(K)
    F = kron(speye(K(i)+1),B);
    theta =@(u) u;
    DAE{i} = set_DAE('E',E, 'A',A, 'F',F, 'theta',theta, 'x0',x0, 'I',I);
end

%% Create high dimensional training set first, project down
P_train = cell(length(K),1);
P_train{end} = rand(K(end)+1,K(end)+200);
for i=length(K)-1:-1:1
    P_train{i} = P_train{i+1}(1:2:end,:);
end
for i=1:length(K)
    P_train{i} = mat2cell(P_train{i}, size(P_train{i},1), ones(size(P_train{i},2),1));
end

%% RB generation through weak greedy algorithm
rb = cell(length(K),1);
for i=1:length(K)
    rb{i} = rb_generation(DAE{i}, P_train{i}, options{i});
end

%% Plot
figure()
semilogy(rb{1}.err_N,'DisplayName',sprintf('K=%i',K(1)))
hold on
for i=2:length(K)
    semilogy(rb{i}.err_N,'DisplayName',sprintf('K=%i',K(i)))
end
title('Maximal greedy training error')
xlabel('RB dimension N')
ylabel('max training error')
legend()
