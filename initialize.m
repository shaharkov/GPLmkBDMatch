function initialize

% reset timer;
tic; dummy = toc;

% set paths
if isempty(whos('global','path_def'))
    fprintf('- Adding toolbox paths\n');
    
    % add common path
    addpath(genpath('code/'));
    addpath(genpath('code_batch/'));
    addpath(genpath('toolbox/'));
    
    if isunix
        % add unix required path
    end
    
    % per-usere path
    user_name = char(java.lang.System.getProperty('user.name'));
    switch user_name
        case {'Shahar Kovalsky', 'shaharko'}
            addpath(genpath('c:/Program Files/Mosek/8/toolbox/r2014a/'));
        case 'tingrangao'
            addpath(genpath('~/Documents/MATLAB/mosek/8/toolbox/r2014a/'));
        case 'trgao10'
            addpath(genpath('~/Documents/MATLAB/mosek/8/tools/platform/osx64x86/bin/'));
        otherwise
            warning('initialize: add your paths here');
    end
    
    global path_def
end
