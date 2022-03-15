function [] = initialise_rb_solver(root)
% RUN THIS FUNCTION RIGHT AT THE START!
% Sets the path for all needed folders.
% Running this function is vital for the rb_solver and its dokumation to
% work.
%
% SYNTAX:
%   initialise_rb_solver
%       if the 'rb_solver' folder is in the current directory.
%   initialise_rb_solver(root)
%       where root is the path to the 'rb_solver' folder.

if ispc
    if nargin == 0
        rb_path = genpath('rb_solver\source');
    else
        root = strrep(root,'/','\');
        rb_path = genpath([root '\rb_solver\source']);
    end
else
    if nargin == 0
        rb_path = genpath('rb_solver/source');
    else
        root = strrep(root,'\','/');
        rb_path = genpath([root '/rb_solver/source']);
    end
end



if isempty(rb_path)
    if nargin == 0
        warning(['Initialisation failed! The rb_solver folder does not exist in ''' pwd '''. ']);
    else
        warning(['Initialisation failed! The rb_solver folder does not exist in ''' root '''. ']);
    end
end

if ispc
    tmp = strsplit(rb_path,';');
else
    tmp = strsplit(rb_path,':');   
end
for i = 2:size(tmp,2)
    if ~contains(path,tmp{i})
        addpath(rb_path);
        fprintf('rb_solver is now ready to use.\n');
        return;
    end
end

fprintf('rb_solver is already inititalised.\n');
    
end