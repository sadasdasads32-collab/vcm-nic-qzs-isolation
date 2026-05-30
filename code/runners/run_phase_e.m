%% run_phase_e.m - Batch run Phase E: Verification & Validation
% Runs all verification scripts sequentially.
% Uses setappdata to persist state across sub-scripts that call `clear;`.

clc;
init_path();

fprintf('===============================================\n');
fprintf('  PHASE E: Verification & Validation Pipeline\n');
fprintf('  Started: %s\n', datestr(now));
fprintf('===============================================\n\n');

root = fileparts(mfilename('fullpath'));
setappdata(0, 'phase_e_root', root);
setappdata(0, 'phase_e_all_passed', true);

% Helper to re-derive directories after possible clear
function [val_dir, stab_dir, opt_dir] = get_dirs()
    root = getappdata(0, 'phase_e_root');
    val_dir = fullfile(root, '..', 'validation');
    stab_dir = fullfile(root, '..', 'stability');
    opt_dir = fullfile(root, '..', 'optimization');
end

function mark_failed()
    setappdata(0, 'phase_e_all_passed', false);
end

%% E1: Harmonic Convergence Verification
[val_dir, stab_dir, opt_dir] = get_dirs();
try
    fprintf('>>> Step 1/4: Verify_Harmonic_Convergence.m\n');
    run(fullfile(val_dir, 'Verify_Harmonic_Convergence.m'));
    fprintf('    DONE.\n\n');
catch ME
    fprintf(2, '    ERROR: %s\n', ME.message);
    fprintf(2, '    %s\n', ME.getReport('basic'));
    mark_failed();
end

%% E2: Bifurcation Classification
[val_dir, stab_dir, opt_dir] = get_dirs();
try
    fprintf('>>> Step 2/4: Run_Bifurcation_Classification.m\n');
    run(fullfile(stab_dir, 'Run_Bifurcation_Classification.m'));
    fprintf('    DONE.\n\n');
catch ME
    fprintf(2, '    ERROR: %s\n', ME.message);
    fprintf(2, '    %s\n', ME.getReport('basic'));
    mark_failed();
end

%% E3: Energy Dissipation & NIC Power Verification
[val_dir, stab_dir, opt_dir] = get_dirs();
try
    fprintf('>>> Step 3/4: Verify_Energy_Dissipation_SCI_v2.m\n');
    run(fullfile(val_dir, 'Verify_Energy_Dissipation_SCI_v2.m'));
    fprintf('    DONE.\n\n');
catch ME
    fprintf(2, '    ERROR: %s\n', ME.message);
    fprintf(2, '    %s\n', ME.getReport('basic'));
    mark_failed();
end

%% E4: Unified Optimization (Data Consistency Check)
[val_dir, stab_dir, opt_dir] = get_dirs();
try
    fprintf('>>> Step 4/4: unified_optimization.m\n');
    fprintf('    This verifies the 62.6%% peak reduction figure.\n');
    run(fullfile(opt_dir, 'unified_optimization.m'));
    fprintf('    DONE.\n\n');
catch ME
    fprintf(2, '    ERROR: %s\n', ME.message);
    fprintf(2, '    %s\n', ME.getReport('basic'));
    mark_failed();
end

all_passed = getappdata(0, 'phase_e_all_passed');
if isempty(all_passed), all_passed = true; end

fprintf('===============================================\n');
if all_passed
    fprintf('  PHASE E COMPLETE - All tests passed!\n');
else
    fprintf('  PHASE E COMPLETE - Some tests had errors (see above).\n');
end
fprintf('  Finished: %s\n', datestr(now));
fprintf('===============================================\n');
