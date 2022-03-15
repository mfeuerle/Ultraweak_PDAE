function [] = log_start(text)
% The logging concept of this framework consist of 5 functions:
%   log_start
%   <a href="matlab: help log_middle">log_middle</a>
%   <a href="matlab: help log_warning">log_warning</a>
%   <a href="matlab: help log_error">log_error</a>
%   <a href="matlab: help log_end">log_end</a>
%
% All functions use the global variable LOG_MODE to control the logging
% level. To change it, just set the value of LOG_MODE to one of these
% values:
%   0:  print all (DEFAULT)
%   1:  print reduced inforamtions
%   2:  print only warnings and errors
%   3:  print only errors
%   4:  print nothing
%
%
% HOW TO USE THE LOGGING FUNCTIONS:
% If you want a function to appear in the logging, just call log_start
% right at the beginning of the function, and <a href="matlab: help
% log_end">log_end</a> at the end of the function. 
% You need both of them!
% The other logging functions can be called between log_start and <a href="matlab: help
% log_end">log_end</a>.
% 
% EXAMPLE:
% function [out] = myFun(in)
%   log_start('starting...');
%   ...
%   log_middle('some additional informations');
%   ...
%   log_end('...finished');
% end

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
end

if isempty(LVL_COUNT)
    LVL_COUNT = 0;
else
    LVL_COUNT = LVL_COUNT + 1;
end

stack = dbstack;

msg = repmat(' | ',[1, LVL_COUNT]);
msg = [msg '<strong>' stack(2).name '</strong>: ' text];

if ~isempty(TIME_RB_START)
    fprintf(['%06.2f ' msg '\n'],toc(TIME_RB_START));
else
    fprintf(['      ' msg '\n']);
end
end

