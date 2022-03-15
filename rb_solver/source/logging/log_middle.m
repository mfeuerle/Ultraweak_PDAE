function [] = log_middle(text)
% log_middle loggs intermediate step with the given message.
%
% For more informations, see <a href="matlab: help log_start">logging</a>.

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

if LOG_MODE > 0
    return;
end

msg = repmat(' | ',[1, LVL_COUNT+1]);
msg = [msg text];

if ~isempty(TIME_RB_START)
    fprintf(['%06.2f ' msg '\n'],toc(TIME_RB_START));
else
    fprintf(['      ' msg '\n']);
end

end

