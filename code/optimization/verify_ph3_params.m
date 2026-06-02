%% verify_ph3_params.m
% Quick验证 Phase 3 最优参数是否真的优于 baseline
% 对比 duibi.m 3000步弧长延拓

clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'lib'));
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'validation'));

% 机械参数
mu=0.2; beta=2; K1=1.0; K2=0.0; U=2.0; Lg=4/9; v=2.5;
alpha1 = v - 2*K1*(1-Lg)/Lg;
alpha2 = beta - 2*K2*(1-Lg)/Lg;
gamma1 = K1/(U^2*Lg^3);
gamma2 = K2/(U^2*Lg^3);

P_base.be1 = 1.0;
P_base.al1 = alpha1 - P_base.be1;
P_base.be2 = alpha2;
P_base.ga1 = gamma1;
P_base.ga2 = gamma2;
P_base.mu  = mu;
P_base.ze1 = 0.05;

global Fw ParamMin ParamMax FixedOmega
Fw = 0.005;
ParamMin = 0.05;
ParamMax = 10.5;
FixedOmega = [];

% ---- Ph3 最优参数 ----
P_opt = P_base;
P_opt.lam   = 0.18;
P_opt.sigma = 1.3728;
P_opt.kap_e = 2.9993;
P_opt.kap_c = 0.6046;

% ---- Baseline (纯机械) ----
P_base2 = P_base;
P_base2.lam   = 0;
P_base2.sigma = 0;
P_base2.kap_e = 0;
P_base2.kap_c = 0;

% sysP = [be1, be2, mu, al1, ga1, ze1, lam, kap_e, kap_c, sigma, ga2]
sysP_opt = [P_opt.be1, P_opt.be2, P_opt.mu, P_opt.al1, P_opt.ga1, P_opt.ze1, ...
            P_opt.lam, P_opt.kap_e, P_opt.kap_c, P_opt.sigma, P_opt.ga2];
sysP_base = [P_base2.be1, P_base2.be2, P_base2.mu, P_base2.al1, P_base2.ga1, P_base2.ze1, ...
             0, 0, 0, 0, P_base2.ga2];

fprintf('=== 验证 Phase 3 最优参数 (3000步弧长) ===\n');
fprintf('sigma=%.4f kap_e=%.4f kap_c=%.4f\n', P_opt.sigma, P_opt.kap_e, P_opt.kap_c);

% --- Ph3 最优参数 FRF ---
tic;
try
    [Om_opt, TF_dB_opt] = arc_length_frf(sysP_opt, 10.0, 'Fw', 0.005, ...
        'Step', -0.01, 'Steps', 3000);
    TF_opt = 10.^(TF_dB_opt/20);
    [TFpk_opt, ipk_opt] = max(TF_opt);
    fprintf('Ph3 最优: TFpk=%.4f (%.1f dB) @ Omega=%.4f\n', ...
        TFpk_opt, 20*log10(TFpk_opt), Om_opt(ipk_opt));
    fprintf('  Om范围: [%.4f, %.4f], Npts=%d\n', min(Om_opt), max(Om_opt), length(Om_opt));
catch e
    fprintf('Ph3 最优: 弧长失败 - %s\n', e.message);
    TFpk_opt = NaN;
end
t1 = toc;

% --- Baseline FRF ---
tic;
try
    [Om_base, TF_dB_base] = arc_length_frf(sysP_base, 10.0, 'Fw', 0.005, ...
        'Step', -0.01, 'Steps', 3000);
    TF_base = 10.^(TF_dB_base/20);
    [TFpk_base, ipk_base] = max(TF_base);
    fprintf('Baseline: TFpk=%.4f (%.1f dB) @ Omega=%.4f\n', ...
        TFpk_base, 20*log10(TFpk_base), Om_base(ipk_base));
    fprintf('  Om范围: [%.4f, %.4f], Npts=%d\n', min(Om_base), max(Om_base), length(Om_base));
catch e
    fprintf('Baseline: 弧长失败 - %s\n', e.message);
    TFpk_base = NaN;
end
t2 = toc;

fprintf('\n耗时: Ph3=%.1fs, Baseline=%.1fs\n', t1, t2);
if ~isnan(TFpk_opt) && ~isnan(TFpk_base)
    reduction = (TFpk_base - TFpk_opt) / TFpk_base * 100;
    fprintf('峰值降低: %.1f%%\n', reduction);
end

% 出对比图
figure('Color','w','Position',[100 100 800 500]);
subplot(1,2,1); hold on; box on; grid on;
set(gca,'XScale','log');
if exist('Om_opt','var')
    semilogx(Om_opt, 20*log10(max(TF_opt,1e-12)), 'b-', 'LineWidth',1.8, ...
        'DisplayName', sprintf('Ph3 opt (TFpk=%.1fdB)', 20*log10(TFpk_opt)));
end
if exist('Om_base','var')
    semilogx(Om_base, 20*log10(max(TF_base,1e-12)), 'r--', 'LineWidth',1.5, ...
        'DisplayName', sprintf('Baseline (TFpk=%.1fdB)', 20*log10(TFpk_base)));
end
yline(0, 'k--'); xlabel('\Omega'); ylabel('T_F (dB)');
title('Ph3 Optimal vs Baseline (3000 steps)');
legend('Location','best');

subplot(1,2,2); hold on; box on; grid on;
if exist('Om_opt','var') && exist('TF_opt','var')
    plot(Om_opt, TF_opt, 'b-', 'LineWidth',1.5);
end
if exist('Om_base','var') && exist('TF_base','var')
    plot(Om_base, TF_base, 'r--', 'LineWidth',1.5);
end
xlabel('\Omega'); ylabel('T_F (linear)');
title('Linear scale');
legend({'Ph3 opt','Baseline'},'Location','best');

fprintf('\nDone.\n');
