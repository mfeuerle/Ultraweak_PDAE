function [] = log_end(supress_or_text)
% log_end indicates the end of a function. Even if you dont want to log any
% message, this method is vital, if you called <a href="matlab: help log_start">log_start</a>
% befor.
%
% INPUT:
%   nothing:    empty line is logged
%         0:    nothing is logged
%    string:    given message is logged
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

if LOG_MODE > 1
    return;
elseif LOG_MODE > 0
    LVL_COUNT = LVL_COUNT - 1;
    return;
else
    if nargin == 0
        msg = repmat(' | ',[1, LVL_COUNT]);
        if ~isempty(TIME_RB_START)
            fprintf(['%06.2f ' msg '\n'],toc(TIME_RB_START));
        else
            fprintf(['      ' msg '\n']);
        end
        %fprintf('\n');
    elseif isnumeric(supress_or_text)
        if supress_or_text ~= 1
            msg = repmat(' | ',[1, LVL_COUNT]);
        if ~isempty(TIME_RB_START)
            fprintf(['%06.2f ' msg '\n'],toc(TIME_RB_START));
        else
            fprintf(['      ' msg '\n']);
        end
        %fprintf('\n');
        end
    else
        msg = [repmat(' | ',[1, LVL_COUNT]) ' ' supress_or_text];
        if ~isempty(TIME_RB_START)
            fprintf(['%06.2f ' msg '\n'],toc(TIME_RB_START));
        else
            fprintf(['      ' msg '\n']);
        end
        %fprintf('\n');
    end
    LVL_COUNT = LVL_COUNT - 1;
    return;
end

end

