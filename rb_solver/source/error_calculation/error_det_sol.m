function L2_norm = error_det_sol(sol, t_analytic,x_analytic, trapz_refinement)

% ----------------------------------------------------------------------
% REFERENCE:
%  E.Beurer, M.Feuerle, N.Reich, K.Urban
%  "An ultraweak variational method for parmeterized linear 
%  differential-algebraic equations"
%  Ulm University, 2022
%  https://doi.org/10.48550/arXiv.2202.12834
%  https://github.com/mfeuerle/Ultraweak_PDAE
% ----------------------------------------------------------------------


if nargin == 4
    f_exact =@(t) interp1(t_analytic,x_analytic,t)';
elseif nargin == 3
    f_exact = t_analytic;
    trapz_refinement = x_analytic;
else
   error('Inputs fehlen'); 
end

x = sol.x;
t = sol.t;

K = size(x,2) -1;

L2_product = 0;
t1 = t(2:end);

switch ndims(x)
    case 2
        x1 = x(:,2:end);
        for k = 1:K
            delta_t = t1(k) - t(k);
            f_approx =@(t_) x(:,k) + (t_-t(k))/delta_t .* (x1(:,k)-x(:,k));
            L2_product = L2_product + integrate_with_trapz(f_exact, f_approx, trapz_refinement, t(k),t1(k));
        end
    case 3
        x1 = x(:,2:end,:);
        for k = 1:K
            delta_t = t1(k) - t(k);
            f_approx =@(t_) x(:,k,2) + (t_-t(k))/delta_t .* (x1(:,k,1)-x(:,k,2));
            L2_product = L2_product + integrate_with_trapz(f_exact, f_approx, trapz_refinement, t(k),t1(k));
        end
    otherwise
        error('was ist hier passier?');
end
L2_norm = sqrt(L2_product);

end

function L2_product = integrate_with_trapz(f_exact, f_approx, tr, t0, t1)
t = linspace(t0,t1,tr+2);
diff = f_exact(t) - f_approx(t);
L2_product = L2_product_discrete(diff,diff,t);
end

function [L2product] = L2_product_discrete(f,g,t)
% calculates the L2 product
%   (f,g)_L2
% of two diskrete functions f,g, given at the times t with the
% trapezoidal rule.

fg = sum(f.*g,1);
L2product = trapz(t,fg);

end

