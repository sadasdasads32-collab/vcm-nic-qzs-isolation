%% verify_peak_reduction.m - Quick verification of peak reduction %
% Computes the exact peak TF for optimized EMSD vs pure mechanical baseline
% at Fw=0.008, to resolve the 61.9% vs 62.6% discrepancy.
clc; init_path();

%% Parameters
mu=0.2; beta=2.0; K1=1.0; K2=0.0; U=2.0; Lg=4/9; v=2.5;
alpha1 = v - 2*K1*(1-Lg)/Lg;
alpha2 = beta - 2*K2*(1-Lg)/Lg;
gamma1 = K1/(U^2*Lg^3);
gamma2 = K2/(U^2*Lg^3);
lam_phys=0.18; theta=sqrt(max(lam_phys,0));
ze1=0.05; be1=1.0; al1=alpha1-be1; be2=alpha2;

sigma_opt=1.150567; kap_e_opt=1.522196; kap_c_opt=0.574336;
Fw_val=0.008;

sysP_opt = [be1,be2,mu,al1,gamma1,ze1,lam_phys,kap_e_opt,kap_c_opt,sigma_opt,gamma2];
sysP_base = [be1,be2,mu,al1,gamma1,ze1,0.0,0.0,0.0,0.0,gamma2];

global FixedOmega Fw
Fw = Fw_val; FixedOmega = [];

Om_vec = logspace(log10(0.2), log10(6.0), 350).';
Nw = length(Om_vec);

fprintf('Peak Reduction Verification\n');
fprintf('Fw=%.4f, 350 freq points in [0.2, 6.0]\n', Fw_val);

%% Sweep for optimized EMSD
fprintf('\nComputing EMSD sweep...\n');
TF_opt = nan(Nw,1);
y_guess = [zeros(15,1); Fw_val];
fail_count = 0;
for j = 1:Nw
    Om = Om_vec(j); FixedOmega = Om;
    try
        y_sol = newton('nondim_temp2', y_guess, sysP_opt);
    catch
        fail_count = fail_count+1;
        if fail_count>12, break; else, continue; end
    end
    xc = y_sol(1:15);
    y_guess = [xc; Fw_val];
    TF_opt(j) = compute_TF_fast(xc, sysP_opt, Om, Fw_val);
end
ok_opt = isfinite(TF_opt) & TF_opt > 0;
peak_opt = max(TF_opt(ok_opt));
peak_opt_dB = 20*log10(peak_opt);
fprintf('EMSD: TF_peak=%.6f (%.2f dB), %d/%d converged\n', peak_opt, peak_opt_dB, nnz(ok_opt), Nw);

%% Sweep for baseline
fprintf('Computing baseline sweep...\n');
TF_base = nan(Nw,1);
y_guess = [zeros(15,1); Fw_val];
fail_count = 0;
for j = 1:Nw
    Om = Om_vec(j); FixedOmega = Om;
    try
        y_sol = newton('nondim_temp2', y_guess, sysP_base);
    catch
        fail_count = fail_count+1;
        if fail_count>12, break; else, continue; end
    end
    xc = y_sol(1:15);
    y_guess = [xc; Fw_val];
    TF_base(j) = compute_TF_fast(xc, sysP_base, Om, Fw_val);
end
ok_base = isfinite(TF_base) & TF_base > 0;
peak_base = max(TF_base(ok_base));
peak_base_dB = 20*log10(peak_base);
fprintf('Baseline: TF_peak=%.6f (%.2f dB), %d/%d converged\n', peak_base, peak_base_dB, nnz(ok_base), Nw);

%% Reduction
reduction = (1 - peak_opt/peak_base) * 100;
fprintf('\n========================================\n');
fprintf('PEAK REDUCTION: %.2f%%\n', reduction);
fprintf('  Baseline TF_peak = %.6f (%.2f dB)\n', peak_base, peak_base_dB);
fprintf('  EMSD TF_peak     = %.6f (%.2f dB)\n', peak_opt, peak_opt_dB);
fprintf('========================================\n');
