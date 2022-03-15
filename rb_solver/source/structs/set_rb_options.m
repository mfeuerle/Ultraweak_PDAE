function options = set_rb_options(varargin)
% set_rb_options settings for the RBM generation.
%
% FIELDS:
%           K: [ positive integer ] time resolution - number of time steps
%         tol: [ positive scalar ]  required tolerance for the RB accuracy
%       N_max: [ positive integer ] max RB dimension
% HomogenMode: [ string ]           Method for homogenization
%                   'sigma0' - substract one hat function in first interval
%                   'const'  - substract constant function on all intervals
%               DEFAULT: sigma0
        
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
    'K          '
    'tol        '
    'Nmax       '
    'HomogenMode'
    'Solver     '
    ];

stack = dbstack;
% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
    help(stack(1).name);
    return;
end

options = struct_handler(Names, varargin);
description =@() help(stack(1).name);
options.help =@() description();

if options.K < 2
    error('K < 2 not allowed.');
end

if isempty(options.HomogenMode)
    options.HomogenMode = 'sigma0';
end

end

