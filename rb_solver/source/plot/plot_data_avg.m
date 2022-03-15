function [t_plot,x_plot] = plot_data_avg(sol)

% ----------------------------------------------------------------------
% REFERENCE:
%  E.Beurer, M.Feuerle, N.Reich, K.Urban
%  "An ultraweak variational method for parmeterized linear 
%  differential-algebraic equations"
%  Ulm University, 2022
%  https://doi.org/10.48550/arXiv.2202.12834
%  https://github.com/mfeuerle/Ultraweak_PDAE
% ----------------------------------------------------------------------

x_plot = 0.5*(sol.x(:,1:end-1,2)+sol.x(:,2:end,1));
t_plot = 0.5*(sol.t(1:end-1)+sol.t(2:end));

end

