function [X_N, R] = GramSchmidt(X_N,R,rb)

% ----------------------------------------------------------------------
% REFERENCE:
%  E.Beurer, M.Feuerle, N.Reich, K.Urban
%  "An ultraweak variational method for parmeterized linear 
%  differential-algebraic equations"
%  Ulm University, 2022
%  https://doi.org/10.48550/arXiv.2202.12834
%  https://github.com/mfeuerle/Ultraweak_PDAE
% ----------------------------------------------------------------------

N = rb.N;
m = rb.DAE.m;
X = rb.pc_data.X;

for i = 1:N
    coef = X.inner_product(X_N(:,N+1),X_N(:,i)) / X.norm(X_N(:,i));
    X_N(:,N+1) = X_N(:,N+1) - coef*X_N(:,i);
    R(:,m+N+2) = R(:,m+N+2) - coef*R(:,m+i+1);
end
norm_N = X.norm(X_N(:,N+1));
X_N(:,N+1) = X_N(:,N+1) / norm_N;
R(:,m+N+2) = R(:,m+N+2) / norm_N;

end