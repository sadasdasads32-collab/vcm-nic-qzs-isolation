%% quick_compare_arc_vs_newton.m
% 快速对比弧长延拓 vs Newton 扫频的 FRF 峰值差异
% 用纯机械基线（lam=0）和一组典型 EMSD 参数分别对比

clc;
init_path();

% --- 机械参数（与 duibi.m 一致）---
mu   = 0.2;
beta = 2.0;
K1   = 1.0;
K2   = 0.0;
U    = 2.0;
Lg   = 4/9;
v    = 2.5;

alpha1 = v    - 2*K1*(1-Lg)/Lg;
alpha2 = beta - 2*K2*(1-Lg)/Lg;
gamma1 = K1/(U^2 * Lg^3);
gamma2 = K2/(U^2 * Lg^3);

be1 = 1.0;
al1 = alpha1 - be1;
be2 = alpha2;

Fw_test = 0.005;

sysP0 = [be1, be2, mu, al1, gamma1, 0.05, 0.18, 1.0, 0.2, 1.0, gamma2];

fprintf('========================================================\n');
fprintf('  弧长延拓 vs Newton 扫频 峰值对比 (Fw=%.4f)\n', Fw_test);
fprintf('========================================================\n\n');

% --- 测试两组参数 ---
test_cases = {
    struct('name', '纯机械基线', 'lam', 0, 'kap_e', 0, 'kap_c', 0, 'sigma', 0), ...
    struct('name', 'EMSD 典型参数', 'lam', 0.18, 'kap_e', 1.5, 'kap_c', 0.6, 'sigma', 1.4)
};

global Fw ParamMin ParamMax
Fw = Fw_test;

for tc = 1:2
    tc_case = test_cases{tc};
    sysP = sysP0;
    sysP(7)  = tc_case.lam;
    sysP(8)  = tc_case.kap_e;
    sysP(9)  = tc_case.kap_c;
    sysP(10) = tc_case.sigma;

    fprintf('--- %s ---\n', tc_case.name);
    fprintf('  参数: lam=%.2f, kap_e=%.2f, kap_c=%.2f, sigma=%.2f\n', ...
        tc_case.lam, tc_case.kap_e, tc_case.kap_c, tc_case.sigma);

    %% 方法 A: 弧长延拓
    ParamMin = 0.05;
    ParamMax = 10.5;

    tic_arc = tic;
    try
        [Om_arc, TF_dB_arc] = arc_length_frf(sysP, 10.0, 'Fw', Fw_test, ...
            'Step', -0.01, 'Steps', 3000);
        TF_lin_arc = 10.^(TF_dB_arc/20);
        peak_arc = max(TF_lin_arc);
        peak_dB_arc = 20*log10(peak_arc);
        [~, idx_pk] = max(TF_lin_arc);
        Om_peak_arc = Om_arc(idx_pk);
        t_arc = toc(tic_arc);
        npts_arc = length(Om_arc);
        fprintf('  弧长延拓: %.1fs, %d 频点, TF_peak=%.4f (%.1f dB) @ Omega=%.4f\n', ...
            t_arc, npts_arc, peak_arc, peak_dB_arc, Om_peak_arc);
    catch e
        fprintf('  弧长延拓: 失败 - %s\n', e.message);
        peak_arc = NaN;
    end

    %% 方法 B: Newton 逐频点扫频
    global FixedOmega
    FixedOmega = [];

    Om_grid = logspace(log10(0.2), log10(6.0), 350).';

    tic_newt = tic;
    TF_newt = nan(size(Om_grid));
    y_guess = [zeros(15,1); Fw_test];
    fail_count = 0;
    for j = 1:length(Om_grid)
        Om = Om_grid(j);
        FixedOmega = Om;
        try
            y_sol = newton('nondim_temp2', y_guess, sysP);
        catch
            fail_count = fail_count + 1;
            if fail_count > 10, break; else, continue; end
        end
        xc = y_sol(1:15);
        y_guess = [xc; Fw_test];
        TF_newt(j) = compute_TF_fast(xc, sysP, Om, Fw_test);
    end
    t_newt = toc(tic_newt);
    ok_newt = ~isnan(TF_newt);
    n_ok = nnz(ok_newt);
    if n_ok > 0
        peak_newt = max(TF_newt(ok_newt));
        peak_dB_newt = 20*log10(peak_newt);
        [~, idx_newt_pk] = max(TF_newt);
        Om_peak_newt = Om_grid(idx_newt_pk);
        fprintf('  Newton扫频: %.1fs, %d/%d 收敛, TF_peak=%.4f (%.1f dB) @ Omega=%.4f\n', ...
            t_newt, n_ok, length(Om_grid), peak_newt, peak_dB_newt, Om_peak_newt);
    else
        fprintf('  Newton扫频: %.1fs, 全部发散\n', t_newt);
        peak_newt = NaN;
    end

    %% 对比
    if ~isnan(peak_arc) && ~isnan(peak_newt)
        ratio = peak_arc / peak_newt;
        dB_diff = peak_dB_arc - peak_dB_newt;
        fprintf('  >>> 峰值比: arc/newt = %.2fx (%.1f dB 差异)\n', ratio, dB_diff);
        if ratio > 1.5
            fprintf('  >>> *** Newton 扫频严重低估真实共振峰！***\n');
        elseif ratio < 0.8
            fprintf('  >>> *** 弧长延拓的峰值低于 Newton 扫频 ***\n');
        else
            fprintf('  >>> 两种方法峰值接近\n');
        end
    end
    fprintf('\n');
end

fprintf('========================================================\n');
fprintf('  结论: 弧长延拓追踪完整解曲线（含共振峰）\n');
fprintf('  Newton 逐频点扫频只追踪低幅值分支，可能错过共振峰\n');
fprintf('========================================================\n');

%% 出对比图
figure('Color','w','Position',[100 100 900 500]);
for tc = 1:2
    subplot(1,2,tc); hold on; box on; grid on;
    tc_case = test_cases{tc};
    sysP = sysP0;
    sysP(7)  = tc_case.lam;
    sysP(8)  = tc_case.kap_e;
    sysP(9)  = tc_case.kap_c;
    sysP(10) = tc_case.sigma;

    ParamMin = 0.05; ParamMax = 10.5;

    try
        [Om_arc, TF_dB_arc] = arc_length_frf(sysP, 10.0, 'Fw', Fw_test, ...
            'Step', -0.01, 'Steps', 3000);
        semilogx(Om_arc, TF_dB_arc, 'b-', 'LineWidth', 1.8, 'DisplayName', 'Arc-length');
    catch
    end

    Om_grid = logspace(log10(0.2), log10(6.0), 350).';
    TF_newt = nan(size(Om_grid));
    FixedOmega = [];
    y_guess = [zeros(15,1); Fw_test];
    for j = 1:length(Om_grid)
        Om = Om_grid(j); FixedOmega = Om;
        try
            y_sol = newton('nondim_temp2', y_guess, sysP);
        catch
            continue;
        end
        xc = y_sol(1:15); y_guess = [xc; Fw_test];
        TF_newt(j) = compute_TF_fast(xc, sysP, Om, Fw_test);
    end
    ok_n = ~isnan(TF_newt);
    TF_dB_newt = 20*log10(max(TF_newt(ok_n), 1e-12));
    semilogx(Om_grid(ok_n), TF_dB_newt, 'r--', 'LineWidth', 1.5, ...
        'DisplayName', 'Newton sweep');

    set(gca, 'XScale', 'log');
    xlabel('\Omega'); ylabel('T_F (dB)');
    title(tc_case.name, 'FontSize', 12);
    yline(0, 'k--', '0 dB');
    legend('Location', 'best');
    xlim([0.1, 10]);
end
sgtitle(sprintf('Arc-length vs Newton Sweep (Fw=%.4f)', Fw_test), 'FontSize', 14);
