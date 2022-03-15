function [t_plot, x_plot] = plot_data(sol)

% ----------------------------------------------------------------------
% REFERENCE:
%  E.Beurer, M.Feuerle, N.Reich, K.Urban
%  "An ultraweak variational method for parmeterized linear 
%  differential-algebraic equations"
%  Ulm University, 2022
%  https://doi.org/10.48550/arXiv.2202.12834
%  https://github.com/mfeuerle/Ultraweak_PDAE
% ----------------------------------------------------------------------

K = size(sol.x,2) - 1;

x_plot = zeros(size(sol.x,1),2*K);
x_plot(:,1:2:2*K-1) = sol.x(:,1:end-1,2);
x_plot(:,2:2:2*K) = sol.x(:,2:end,1);

t_plot = [sol.t(1); repelem(sol.t(2:end-1),2); sol.t(end)];
end

