function init_path()
%% init_path - Add all project subdirectories to MATLAB path
% Run this once at the start of each MATLAB session, or call
% init_path() at the top of any run script.

    root = fileparts(mfilename('fullpath'));
    if isempty(root)
        root = pwd;
    end

    folders = {'lib', 'frf_sweep', 'stability', 'boa', ...
               'operator', 'optimization', 'validation', 'figures', 'runners'};

    fprintf('Adding project paths...\n');
    for i = 1:numel(folders)
        p = fullfile(root, folders{i});
        if exist(p, 'dir')
            addpath(p);
            fprintf('  + %s\n', folders{i});
        end
    end
    fprintf('Done. All paths added.\n');
end
