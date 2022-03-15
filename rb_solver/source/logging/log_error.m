function [] = log_error(text)
% Not ready yet, dont use it

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
global LOG_MODE;

if LOG_MODE > 3
    return;
end

stack = dbstack;

msg = repmat(' | ',[1, LVL_COUNT]);
msg = [msg '<strong>' stack(2).name '</strong>: ' text];

if ~isempty(TIME_RB_START)
    error(['%06.2f ' msg '\n'],toc(TIME_RB_START));
else
    error(['      ' msg '\n']);
end
end

