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
% for local (fast) computation
K  = [5,10,15,20];
% as in the paper
% K = [100,150,200,250,300];

%DAE
nx = 15;
ny = 15;
m = 1;

%RB
% fast computation
resolutions = 2:2:18;
% paper
% resolutions = [10 30 50 70 90 130 170 210];
reduced_dimensions = cell(length(K),1);
for i=1:length(K)
    reduced_dimensions{i} = resolutions(resolutions<K(i));
end

% Validation with Sine
sine_periods = 1:15;

% Validation with Step function
step_length = 0.2;
step_starts = 0:0.1:0.8;

%% RBM Options
HomogenMode = 'sigma0'; % Way how the homogenization of x0 is applied
solver  = set_solver('direct');
options = cell(length(K),1);
for i=1:length(K)
    options{i} = set_rb_options('K',           K(i),...
                                'HomogenMode', HomogenMode,...
                                'Solver',      solver);
end
                     
%% Time Discretisation
I = [0 1]; % interval
t = cell(length(K),1);
for i=1:length(K)
    t{i} = linspace(I(1),I(2), K(i)+1);
end

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
theta =@(u) u;
for i=1:length(K)
    F = kron(speye(K(i)+1),B);
    DAE{i} = set_DAE('E',E, 'A',A, 'F',F, 'theta',theta, 'x0',x0, 'I',I);
end

%% Precalculation
pc_data = cell(length(K),1);
rb_init = cell(length(K),1);
for i=1:length(K)
    pc_data{i} = precalculate_data(DAE{i}, options{i});
    rb_init{i} = initialising_rb(DAE{i}, pc_data{i}, [], options{i});
end

%% RBM generation for different time resolutions of the control (Sine)
p_sine = length(sine_periods);
p_step = length(step_starts);
p_all = p_sine + p_step;
P_validation = cell(length(K),1);
for j=1:length(K)
    P_validation{j} = zeros(K(j)+1,p_all);
    for i=1:p_all
        if i <= p_sine
            P_validation{j}(:,i) = sin(sine_periods(i)*2*pi*t{j}/(I(2)-I(1)));
        else
            step_start = step_starts(i-p_sine);
            P_validation{j}(:,i) = 0.5*((t{j}>=step_start)&(t{j}<(step_start+step_length))) + 0.5;
        end
    end
end

err= cell(length(K),1);
err_rel = cell(length(K),1);

for i=1:length(K)
    err{i} = zeros(length(reduced_dimensions{i}),p_all);
    err_rel{i} = zeros(length(reduced_dimensions{i}),p_all);
end

sol = cell(length(K),p_all);
norm_sol = zeros(length(K),p_all);
for i=1:length(K)
    for j = 1:p_all
        sol{i,j} = solver_DAE(DAE{i}, P_validation{i}(:,j), options{i}, pc_data{i});
        norm_sol(i,j) = pc_data{i}.X.norm(sol{i,j}.x_vec);
    end
end

rbs = cell(length(K),length(reduced_dimensions{end}));
for j=1:length(K)
    for i = 1:length(reduced_dimensions{j})
        P_train = discrete_piecewise_linear_basis(reduced_dimensions{j}(i),t{j});
        rb_init_temp = rb_init{j};
        rb_init_temp.P_train = P_train;

        rb = rb_generation_simple(DAE{j}, P_train, options{j}, rb_init_temp);
        rbs{j,i} = rb;

        err_temp = zeros(1,p_all);
        err_rel_temp = zeros(1,p_all);
        for k = 1:p_all
            sol_rb = solver_rb(rb, P_validation{j}(:,k));
            err_temp(k) = error_rb(sol_rb, sol{j,k});
            err_rel_temp(k) = err_temp(k) / norm_sol(j,k);
        end
        err{j}(i,:) = err_temp';
        err_rel{j}(i,:) = err_rel_temp';
    end
end

%% Plotting
figure()
for j = 1:length(K)
    values = max(err_rel{j},[],2);
    loglog(resolutions(1:length(values)),values,'s-',...
        'DisplayName',['K = ' num2str(K(j))])
    hold on
end
xlabel('reduced control dimension K_u')
ylabel('max relative error')
title('error for validation set with sine and discontinuous functions')
legend()
