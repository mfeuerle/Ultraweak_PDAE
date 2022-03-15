function problem = set_DAE(varargin)
% set_DAE creates a linear Differential-algebraic system with parameter
% dependen right hand side with affine decomposition. Given a parameter mu,
% we seek for the solution x : I -> R^n
%
%      E*x'(t;mu) - A*x(t;mu) = F(t)*theta(mu)      for all t in I
%      x(0;mu)  = x0
%
% FIELDS:
%     E: [ matrix R^n*n ]           probably singular, but pencil {E,A} regular
%     A: [ matrix R^n*n ]           probably singular, but pencil {E,A} regular
%     F: [ func, I -> R^n*m ]       affine decomposition of the rhs
%        [ matrix R^(n*(K+1))*m ]   function handle or discrete
% theta: [ func, ParamSpace -> R^m ]affine decomposition of the rhs
%    x0: [ vector R^n ]             initial state
%     I: [ positive scalar tupel ]  time interval
%     n: [ positive integer ]       state dimension          | automatically set
%     m: [ positive integer ]       affine decomp. dimension | automatically set

% ----------------------------------------------------------------------
% REFERENCE:
%  E.Beurer, M.Feuerle, N.Reich, K.Urban
%  "An ultraweak variational method for parmeterized linear 
%  differential-algebraic equations"
%  Ulm University, 2022
%  https://doi.org/10.48550/arXiv.2202.12834
%  https://github.com/mfeuerle/Ultraweak_PDAE
% ----------------------------------------------------------------------

Names = [
    'E    '
    'A    '
    'F    '
    'theta'
    'x0   '
    'I    '
    'n    '
    'm    '
    ];

stack = dbstack;
% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
    help(stack(1).name);
    return;
end

problem = struct_handler(Names, varargin);
description =@() help(stack(1).name);
problem.help =@() description();

if 2 ~= length(problem.I) || problem.I(2) <= problem.I(1)
    error('I must be a time interval');
end

if isa(problem.F,'function_handle')
    [n,m] = size(problem.F(problem.I(1)));
else
    m = size(problem.F,2);
    n = size(problem.A,1);
end

if n ~= size(problem.E,1) || n ~= size(problem.E,2)
   error('E is not in R^n*n');
end
if n ~= size(problem.A,1) || n ~= size(problem.A,2)
   error('A is not in R^n*n');
end

x0 = problem.x0;
if ~(n == size(x0,1) && 1 == size(x0,2))
    if 1 == size(x0,1) && n == size(x0,2)
        problem.x0 = x0';
    else
        error('x0 is not in R^n');
    end
end

problem.n = n;
problem.m = m;

end

