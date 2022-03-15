function [] = log_warning(text)
% log_warning loggs a warning, with the given message.
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

if LOG_MODE > 2
    return;
end

stack = dbstack;

%msg1 = repmat(' | ',[1, LVL_COUNT]);
%msg2 = ['<strong>' stack(2).name '</strong>: ' text];
msg1 = repmat(' | ',[1, LVL_COUNT+1]);
msg2 = ['<strong>' text '</strong>'];
msg2 = text;

if ~isempty(TIME_RB_START)
    %fprintf(['[\b%06.2f ' msg1 msg2 '\n]\b'],toc(TIME_RB_START));
    fprintf(['[\b%06.2f]\b ' msg1 '[\b' msg2 '\n]\b'],toc(TIME_RB_START));
else
    fprintf(['[\b      ' msg2 '\n]\b']);
end
end

