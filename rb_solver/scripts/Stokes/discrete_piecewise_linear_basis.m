function [U] = discrete_piecewise_linear_basis(base_dim,t)
% Calculate a piecewise linear basis via hat functions for a space
% dimension and time grid

% ----------------------------------------------------------------------
% REFERENCE:
%  E.Beurer, M.Feuerle, N.Reich, K.Urban
%  "An ultraweak variational method for parmeterized linear 
%  differential-algebraic equations"
%  Ulm University, 2022
%  https://doi.org/10.48550/arXiv.2202.12834
%  https://github.com/mfeuerle/Ultraweak_PDAE
% ----------------------------------------------------------------------

if base_dim == 1
   U = {zeros(length(t),1), ones(length(t),1)};
   return
end

points = linspace(t(1),t(end),base_dim);

delta_t = (t(end) - t(1)) / (base_dim-1);

index = ( t' >= points(1:end-1) ) & ( t' < points(2:end) );

U = [index .* (points(2:end)-t'), zeros(length(t),1)]...
  + [zeros(length(t),1), index .* (t'-points(1:end-1))];
U = U / delta_t;
U(end,end) = 1;

U = [zeros(length(t),1), U];

U = mat2cell(U,length(t),ones(base_dim+1,1));
end

