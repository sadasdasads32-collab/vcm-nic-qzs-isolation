%% run_full_optimization.m - Wrapper to run unified optimization
% Ensures correct path setup before running the optimization

clc;
init_path();

diary_file = fullfile(fileparts(mfilename('fullpath')), '..', 'logs', 'opt_results.log');
diary(diary_file);

fprintf('========================================\n');
fprintf('  Full Unified Optimization Pipeline\n');
fprintf('  Started: %s\n', datestr(now));
fprintf('========================================\n\n');

% Run the optimization
unified_optimization;

diary off;
fprintf('\n========================================\n');
fprintf('  Optimization complete!\n');
if exist('diary_file', 'var')
    fprintf('  Results saved to: %s\n', diary_file);
else
    fprintf('  Results saved to: data/ & logs/\n');
end
fprintf('========================================\n');
