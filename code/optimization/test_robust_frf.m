%% test_robust_frf.m — Quick test of arc_length_frf_robust vs Ph3 params
clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'lib'));

mu=0.2; beta=2; K1=1.0; K2=0.0; U=2.0; Lg=4/9; v=2.5;
alpha1=v-2*K1*(1-Lg)/Lg; alpha2=beta-2*K2*(1-Lg)/Lg;
gamma1=K1/(U^2*Lg^3); gamma2=K2/(U^2*Lg^3);
be1=1.0; al1=alpha1-be1; be2=alpha2;

global Fw ParamMin ParamMax FixedOmega
Fw=0.005; FixedOmega=[]; ParamMin=0.05; ParamMax=10.5;

% Ph3 Ooptimal params
sysP_opt = [be1, be2, mu, al1, gamma1, 0.05, 0.18, 2.9993, 0.6046, 1.3728, gamma2];

% Baseline
sysP_base = [be1, be2, mu, al1, gamma1, 0.05, 0, 0, 0, 0, gamma2];

fprintf('=== arc_length_frf_robust (Budget=3000) ===\n');
tic;
[Om_opt, TF_dB_opt] = arc_length_frf_robust(sysP_opt, 10.0, 'Fw', 0.005, 'Budget', 3000, 'Verbose', true);
TF_opt = 10.^(TF_dB_opt/20);
[TFpk_opt, ipk] = max(TF_opt);
fprintf('Ph3: TFpk=%.4f (%.1f dB) @ Om=%.4f, OmRange=[%.4f,%.4f], Npts=%d\n', ...
    TFpk_opt, 20*log10(TFpk_opt), Om_opt(ipk), min(Om_opt), max(Om_opt), length(Om_opt));
t1=toc;

fprintf('\n=== arc_length_frf_robust (Budget=3000) Baseline ===\n');
tic;
[Om_base, TF_dB_base] = arc_length_frf_robust(sysP_base, 10.0, 'Fw', 0.005, 'Budget', 3000, 'Verbose', true);
TF_base = 10.^(TF_dB_base/20);
[TFpk_base, ipk2] = max(TF_base);
fprintf('Base: TFpk=%.4f (%.1f dB) @ Om=%.4f, OmRange=[%.4f,%.4f], Npts=%d\n', ...
    TFpk_base, 20*log10(TFpk_base), Om_base(ipk2), min(Om_base), max(Om_base), length(Om_base));
t2=toc;

if ~isempty(Om_opt) && ~isempty(Om_base)
    reduction = (TFpk_base - TFpk_opt)/TFpk_base*100;
    fprintf('\nReduction: %.1f%% (%.1fs+%.1fs)\n', reduction, t1, t2);
end

% Quick plot
figure('Color','w'); hold on; box on; grid on;
if ~isempty(Om_opt)
    semilogx(Om_opt, 20*log10(max(TF_opt,1e-12)), 'b-', 'LineWidth',1.5, ...
        'DisplayName', sprintf('Ph3 opt (%.1f dB)', 20*log10(TFpk_opt)));
end
if ~isempty(Om_base)
    semilogx(Om_base, 20*log10(max(TF_base,1e-12)), 'r--', 'LineWidth',1.5, ...
        'DisplayName', sprintf('Baseline (%.1f dB)', 20*log10(TFpk_base)));
end
yline(0, 'k--'); set(gca,'XScale','log');
xlabel('\Omega'); ylabel('T_F (dB)');
title('arc\_length\_frf\_robust: Ph3 vs Baseline');
legend('Location','best');
fprintf('Done.\n');
