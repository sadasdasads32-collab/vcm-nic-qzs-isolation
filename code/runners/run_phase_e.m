%% run_phase_e.m - Batch run Phase E: Verification & Validation
% Runs all verification scripts sequentially

clc;
init_path();

fprintf('===============================================\n');
fprintf('  PHASE E: Verification & Validation Pipeline\n');
fprintf('  Started: %s\n', datestr(now));
fprintf('===============================================\n\n');

root = fileparts(mfilename('fullpath'));
val_dir = fullfile(root, '..', 'validation');
stab_dir = fullfile(root, '..', 'stability');
opt_dir = fullfile(root, '..', 'optimization');

all_passed = true;

%% E1: Harmonic Convergence Verification
try
    fprintf('>>> Step 1/4: Verify_Harmonic_Convergence.m\n');
    run(fullfile(val_dir, 'Verify_Harmonic_Convergence.m'));
    fprintf('    DONE.\n\n');
catch ME
    fprintf(2, '    ERROR: %s\n', ME.message);
    fprintf(2, '    %s\n', ME.getReport('basic'));
    all_passed = false;
end

%% E2: Bifurcation Classification
try
    fprintf('>>> Step 2/4: Run_Bifurcation_Classification.m\n');
    run(fullfile(stab_dir, 'Run_Bifurcation_Classification.m'));
    fprintf('    DONE.\n\n');
catch ME
    fprintf(2, '    ERROR: %s\n', ME.message);
    fprintf(2, '    %s\n', ME.getReport('basic'));
    all_passed = false;
end

%% E3: Energy Dissipation & NIC Power Verification
try
    fprintf('>>> Step 3/4: Verify_Energy_Dissipation_SCI_v2.m\n');
    run(fullfile(val_dir, 'Verify_Energy_Dissipation_SCI_v2.m'));
    fprintf('    DONE.\n\n');
catch ME
    fprintf(2, '    ERROR: %s\n', ME.message);
    fprintf(2, '    %s\n', ME.getReport('basic'));
    all_passed = false;
end

%% E4: Unified Optimization (Data Consistency Check)
try
    fprintf('>>> Step 4/4: unified_optimization.m\n');
    fprintf('    This verifies the 62.6%% peak reduction figure.\n');
    run(fullfile(opt_dir, 'unified_optimization.m'));
    fprintf('    DONE.\n\n');
catch ME
    fprintf(2, '    ERROR: %s\n', ME.message);
    fprintf(2, '    %s\n', ME.getReport('basic'));
    all_passed = false;
end

fprintf('===============================================\n');
if all_passed
    fprintf('  PHASE E COMPLETE - All tests passed!\n');
else
    fprintf('  PHASE E COMPLETE - Some tests had errors (see above).\n');
end
fprintf('  Finished: %s\n', datestr(now));
fprintf('===============================================\n');
