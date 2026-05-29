%% run_phase_d.m - Batch run Phase D: Figure Regeneration
% Runs all three figure generation scripts sequentially

clc;
init_path();

fprintf('===============================================\n');
fprintf('  PHASE D: Figure Regeneration Pipeline\n');
fprintf('  Started: %s\n', datestr(now));
fprintf('===============================================\n\n');

try
    fprintf('>>> Step 1/3: Generate_All_Journal_Figures.m\n');
    fprintf('    This generates the core journal figures (EPS format)\n');
    run(fullfile(fileparts(mfilename('fullpath')), '..', 'figures', 'Generate_All_Journal_Figures.m'));
    fprintf('    DONE.\n\n');
catch ME
    fprintf(2, '    ERROR: %s\n', ME.message);
end

try
    fprintf('>>> Step 2/3: Generate_Chapter4_Figures.m\n');
    fprintf('    This generates Chapter 4 K(Omega) figures (PDF format)\n');
    run(fullfile(fileparts(mfilename('fullpath')), '..', 'figures', 'Generate_Chapter4_Figures.m'));
    fprintf('    DONE.\n\n');
catch ME
    fprintf(2, '    ERROR: %s\n', ME.message);
end

try
    fprintf('>>> Step 3/3: Generate_All_Chapter5_6_Figures.m\n');
    fprintf('    This generates Chapters 5-6 FRF/Stability figures (PDF format)\n');
    run(fullfile(fileparts(mfilename('fullpath')), '..', 'figures', 'Generate_All_Chapter5_6_Figures.m'));
    fprintf('    DONE.\n\n');
catch ME
    fprintf(2, '    ERROR: %s\n', ME.message);
end

fprintf('===============================================\n');
fprintf('  PHASE D COMPLETE: %s\n', datestr(now));
fprintf('===============================================\n');

% List generated files
out_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'output', 'journal_figures');
if exist(out_dir, 'dir')
    files = dir(fullfile(out_dir, '*'));
    fprintf('\nAll files in %s:\n', out_dir);
    for k = 1:length(files)
        if ~files(k).isdir
            fprintf('  %s  (%.1f KB)\n', files(k).name, files(k).bytes/1024);
        end
    end
end
