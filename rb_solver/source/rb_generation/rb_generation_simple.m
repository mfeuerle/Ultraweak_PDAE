function [rb] = rb_generation_simple(DAE, P_train, options, rb)
% Simple RB generation, just adding all parameters to the basis.

% ----------------------------------------------------------------------
% REFERENCE:
%  E.Beurer, M.Feuerle, N.Reich, K.Urban
%  "An ultraweak variational method for parmeterized linear 
%  differential-algebraic equations"
%  Ulm University, 2022
%  https://doi.org/10.48550/arXiv.2202.12834
%  https://github.com/mfeuerle/Ultraweak_PDAE
% ----------------------------------------------------------------------

global TIME_RB_START;
global LVL_COUNT;
TIME_RB_START  = tic;
LVL_COUNT = [];
warning ('on','all');

log_start(['constructing a reduced basis with N = ' num2str(length(P_train)) '...']);

%% Pre-calculating and initialising
log_middle('initialising...');
if nargin < 4 || isempty(rb)
    pc_data = precalculate_data(DAE, options);
    rb = initialising_rb(DAE, pc_data, P_train, options);
end

%% Full Iteration
N = length(P_train);
log_middle('starting rb generation...');
for i = 1:N
    log_middle(sprintf('%04.1f%%%%, %3u of %u snapshots done.',(i-1)/N*100,i-1,N));
    rb = add_new_basis(rb,i);
end
        
log_end('100.0%%,rb generation finished.');
end

