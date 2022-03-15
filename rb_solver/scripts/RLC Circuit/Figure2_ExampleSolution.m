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

for smooth = [0, 1]
    %% RLC Circuit Problem
    L = 1;  % inductivity in Henry
    C = 1;  % capacity in Farad
    R = 1;  % resistance in Ohm

    E = [1   0   0   0;
         0   1   0   0;
         0   0   0   0;
         0   0   0   0];
    A = [0   0 1/L   0;
         1/C 0   0   0;
         R   0   0  -1;
         0   1   1   1];
    B = [0;  0;  0; -1];

    %% Voltage source and corresponding analytical solution (for smooth case)
    U_hat = 5; % voltage amplitude in Volt
    phi_U = 0; % phase shift at t0 in rad
    if smooth
        omega = 1/sqrt(L*C); % frequency in 1/s
        f_Vs =@(t) U_hat*cos(omega*t);
    else
        omega = 120;         % frequency in 1/s
        f_Vs =@(t) U_hat*sign(cos(omega*t));
    end

    [x_analy,~,dx_analy] = get_analytic_solution(R,L,C,U_hat,phi_U,omega);

    %% Time intervall and discretisation
    t0 = 0; % starting point
    T  = t0 + 2*2*pi/omega; % end
    I = [t0 T]; % interval

    K  = 100;
    t = linspace(I(1),I(2), K+1);
    
    %% Assamble affine decomposition of the rhs
    F = kron(speye(K+1),B);
    theta =@(u) u(t');
    
    %% Assamble the full problem
    x0 = x_analy(I(1));
    DAE = set_DAE('E',E, 'A',A, 'F',F, 'theta',theta, 'x0',x0, 'I',I);

    %% Calculate Petrov-Galerking approximation
    HomogenMode = 'const'; % Way how the homogenization of x0 is applied
    solver  = set_solver('direct'); % solver for the LGS
    options = set_rb_options('K',K, 'HomogenMode',HomogenMode, 'Solver',solver);

    sol = solver_DAE(DAE, f_Vs, options);
    
    %% Calculate ode15i solution
    f15i =@(t,x,dx) E*dx - A*x - B*f_Vs(t);
    x0 = x_analy(I(1));
    dx0 = dx_analy(I(1));
    [t_ode, x_ode] = ode15i(f15i,I,x0,dx0);

    %% Plot 
    [t_pg,x_pg] = plot_data_avg(sol);
    
    xLabels = {'current [A]', 'voltage at capacitor [V]',...
               'voltage at inductor [V]', 'voltage at resistor [V]'};
    figure()
    for i = 1:4
        subplot(2,2,i)
        hold on 
        if smooth
            x = x_analy(t);
            plot(t,x(i,:) , 'k', 'DisplayName', 'analytical solution')
        end
        plot(t_pg, x_pg(i,:),'r', 'DisplayName', 'ultraweak')
        plot(t_ode, x_ode(:,i), 'b--', 'DisplayName', 'ode15i')

        xlabel(xLabels{i})
        xlim([I(1), I(2)])
        legend()
        
    end
    if smooth
        sgtitle('smooth right hand side')
    else
        sgtitle('discontinuous right hand side')
    end
end