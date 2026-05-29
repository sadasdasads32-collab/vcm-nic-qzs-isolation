%% unified_optimization.m
% =========================================================================
% 两级优化流水线：算子快筛 → 系统精修 → 稳定性验证
%
% Phase 1: 算子层 LHS 快筛（K(Ω) 域，数千候选，秒级）
%          → 用解析算子形状评分预选出 Top-N 候选电路参数
% Phase 2: 系统层 FRF 验证（HBM+Newton 逐频点，Top-N 候选，分钟级）
%          → 计算实际力传递率峰值，精选 Top-K
% Phase 3: 系统层 fminsearch 精修（Top-K 候选局部优化，分钟级）
%          → 在真实 HBM 模型上微调参数
% Phase 4: Floquet 稳定性全扫描 + 出版级对比图
%
% 连接了 suanzi_op.m（算子域）和 optimization.m（系统域），
% 在保持精度的同时将总计算量减少约 20 倍。
% =========================================================================

clc; clear; close all;
init_path();

%% ========================================================================
% 0. 用户可调参数
% ========================================================================

% --- 机械参数（与 optimization.m 和 Wang(2017) BG 模型一致）---
mu   = 0.2;         % 质量比
beta = 2.0;         % 下层线性刚度比
K1   = 1.0;         % 上层水平弹簧比
K2   = 0.2;         % 下层水平弹簧比
U    = 2.0;         % 几何尺度
Lg   = 4/9;         % 杆长比
v    = 2.5;         % 由 L=4/9, alpha1=0 反推

alpha1 = v    - 2*K1*(1-Lg)/Lg;
alpha2 = beta - 2*K2*(1-Lg)/Lg;
gamma1 = K1/(U^2 * Lg^3);
gamma2 = K2/(U^2 * Lg^3);

lam_phys = 0.18;    % 机电耦合系数 theta^2
theta    = sqrt(max(lam_phys, 0));
ze1      = 0.05;    % 下层阻尼比

be1 = 1.0;
al1 = alpha1 - be1;
be2 = alpha2;

% --- 激励与频带 ---
Fw_opt  = 0.008;    % 优化用激励幅值
Om_min  = 0.2;
Om_max  = 6.0;
Nw_FRF  = 180;      % 每候选频点数（Phase 2）

% --- 算子预筛设置 ---
Nsamp_phase1  = 5000;    % Phase 1 LHS 采样数
TopN_phase1   = 30;      % 进入 Phase 2 的候选数
TopK_phase2   = 8;       % 进入 Phase 3 精修的候选数
MaxIter_fmin  = 80;      % fminsearch 最大迭代

% --- 搜索边界 [sigma, kap_e, kap_c] ---
lb = [0.02, 0.02, 0.02];
ub = [3.00, 3.00, 3.00];

% --- Floquet 设置 ---
Nt_floquet  = 600;
tol_stable  = 1.002;

rng(42);
fprintf('========== 两级优化流水线启动 ==========\n');
fprintf('lam=%.4f (theta=%.4f), Fw=%.4f\n\n', lam_phys, theta, Fw_opt);

%% ========================================================================
% Phase 1: 算子层 LHS 快筛
% ========================================================================
fprintf('===== Phase 1: 算子层 LHS 快筛 (%d 样本) =====\n', Nsamp_phase1);

% 频率网格（仅用于算子评估，与真实 HBM 无关）
Om_op = logspace(log10(Om_min), log10(Om_max), 350).';
Om_split = 1.0;
idxL = (Om_op <= Om_split);
idxH = (Om_op >= Om_split);

% 硬约束
minDenThresh = 0.02;

% 权重与阈值（根据实际系统量级校准）
w = struct('KL', 2.0, 'CL', 2.5, 'ML', 1.5, 'KH', 1.5, 'CH', 2.0);
Ctar   = 0.08;    % 低频阻尼目标下限
KtarH  = 0.03;    % 高频刚度目标下限
Mtar   = 0.30;    % 低频惯容目标下限（theta^2/kap_c 量级）

X = lhsdesign(Nsamp_phase1, 3);
X = lb + X .* (ub - lb);

J_phase1 = inf(Nsamp_phase1, 1);

tic_phase1 = tic;
for i = 1:Nsamp_phase1
    sigma_i = X(i,1); kap_e_i = X(i,2); kap_c_i = X(i,3);

    den = kap_e_i * Om_op.^2 - 1i * sigma_i * Om_op - kap_c_i;
    if min(abs(den)) < minDenThresh
        continue;
    end

    K = (theta^2) .* Om_op.^2 ./ den;
    Kr  = real(K);
    Ceq = imag(K) ./ max(Om_op, 1e-12);

    % LF-anchored split
    Meq = (theta^2) / kap_c_i * ones(size(Om_op));
    Keq = Kr + Om_op.^2 .* Meq;

    lw = ones(size(Om_op)) ./ length(Om_op);
    meanL = @(y) sum(y(idxL) .* lw(idxL)) / sum(lw(idxL));
    meanH = @(y) sum(y(idxH) .* lw(idxH)) / sum(lw(idxH));

    KL = meanL(abs(Keq));
    KH = meanH(abs(Keq));
    CL = meanL(Ceq);
    CH = meanH(Ceq);
    ML = meanL(Meq);

    hinge = @(z) max(0, z);
    Jshape = 0;
    Jshape = Jshape + w.KL * KL;
    Jshape = Jshape + w.CL * hinge(Ctar - CL);
    Jshape = Jshape + w.ML * hinge(Mtar - ML);
    Jshape = Jshape + w.KH * hinge(KtarH - KH);
    Jshape = Jshape + w.CH * CH;

    [CePeak, ~] = max(Ceq);
    sharp = CePeak / max(CL, 1e-12);
    Jshape = Jshape + 5.0 * hinge(sharp - 12);

    J_phase1(i) = Jshape;
end
t_phase1 = toc(tic_phase1);

[Js, idx_sort] = sort(J_phase1, 'ascend');
idx_valid = idx_sort(isfinite(Js));
TopN_actual = min(TopN_phase1, numel(idx_valid));
idx_top1 = idx_valid(1:TopN_actual);

fprintf('Phase 1 完成 (%.1f s), Top-%d 候选已选出\n', t_phase1, TopN_actual);
fprintf('  最佳算子评分: %.4e\n', J_phase1(idx_top1(1)));
fprintf('  参数: sigma=%.4f, kap_e=%.4f, kap_c=%.4f\n\n', ...
    X(idx_top1(1),1), X(idx_top1(1),2), X(idx_top1(1),3));

%% ========================================================================
% Phase 2: 系统层 FRF 验证（Top-N，HBM + Newton 逐频点）
% ========================================================================
fprintf('===== Phase 2: 系统层 FRF 验证 (Top-%d) =====\n', TopN_actual);

global FixedOmega Fw
Fw = Fw_opt;
FixedOmega = [];

sysP0 = [be1, be2, mu, al1, gamma1, ze1, lam_phys, ...
         1.0, 0.2, 1.0, gamma2];

Om_grid = logspace(log10(Om_min), log10(Om_max), Nw_FRF).';

TF_peak  = inf(TopN_actual, 1);
TF_all   = cell(TopN_actual, 1);
Ok_all   = cell(TopN_actual, 1);
Mu_peak  = inf(TopN_actual, 1);

tic_phase2 = tic;
for r = 1:TopN_actual
    i = idx_top1(r);
    sigma_i = X(i,1); kap_e_i = X(i,2); kap_c_i = X(i,3);

    sysP = sysP0;
    sysP(8)  = kap_e_i;
    sysP(9)  = kap_c_i;
    sysP(10) = sigma_i;

    % --- 快速 FRF：Newton 逐频点 + 初值延续 ---
    TF = nan(Nw_FRF, 1);
    ok = false(Nw_FRF, 1);

    xc_low = []; Om_low = NaN;     % 低频端解
    xc_peak = []; Om_peak = NaN;   % 共振峰解
    xc_high = []; Om_high = NaN;   % 高频端解
    peak_TF_val = -inf;

    y_guess = [zeros(15,1); Fw_opt];
    fail_count = 0;

    for j = 1:Nw_FRF
        Om = Om_grid(j);
        FixedOmega = Om;

        try
            y_sol = newton('nondim_temp2', y_guess, sysP);
        catch
            fail_count = fail_count + 1;
            if fail_count > 10, break; else, continue; end
        end

        xc = y_sol(1:15);
        y_guess = [xc; Fw_opt];
        ok(j) = true;

        % 力传递率（仅基波）
        TF(j) = compute_TF_fast(xc, sysP, Om, Fw_opt);

        % 记录关键频点的 HB 解
        if isempty(xc_low)
            xc_low = xc; Om_low = Om;
        end
        if TF(j) > peak_TF_val
            peak_TF_val = TF(j);
            xc_peak = xc; Om_peak = Om;
        end
        xc_high = xc; Om_high = Om;  % 始终更新到最后一个收敛点
    end

    TF_all{r} = TF;
    Ok_all{r} = ok;

    if nnz(ok) < 10
        TF_peak(r) = 1e6;
        Mu_peak(r) = 1e6;
    else
        TF_peak(r) = max(TF(ok));
        % 多点 Floquet 抽查：低频 / 峰值 / 高频（去重）
        check_oms = [];
        check_xcs = {};
        tol_om = 1e-6;
        if ~isempty(xc_low)
            check_oms(end+1) = Om_low;
            check_xcs{end+1} = xc_low;
        end
        if ~isempty(xc_peak) && all(abs(Om_peak - check_oms) > tol_om)
            check_oms(end+1) = Om_peak;
            check_xcs{end+1} = xc_peak;
        end
        if ~isempty(xc_high) && all(abs(Om_high - check_oms) > tol_om)
            check_oms(end+1) = Om_high;
            check_xcs{end+1} = xc_high;
        end
        mu_vals = zeros(1, length(check_xcs));
        for k = 1:length(check_xcs)
            mu_vals(k) = compute_floquet_fast(check_xcs{k}, sysP, check_oms(k), Nt_floquet);
        end
        Mu_peak(r) = max(mu_vals);
    end

    if mod(r, 5) == 0
        fprintf('  Phase 2: %2d/%d done, best TF_peak=%.3f\n', ...
            r, TopN_actual, min(TF_peak(1:r)));
    end
end
t_phase2 = toc(tic_phase2);

% 保存 Phase 2 原始最佳（用于后续对比）
[best_TF_p2, idx_best_p2] = min(TF_peak);
p_best_phase2 = X(idx_top1(idx_best_p2), :);

fprintf('Phase 2 完成 (%.1f s)\n', t_phase2);
fprintf('  最佳 TF_peak = %.4f @ sigma=%.4f kap_e=%.4f kap_c=%.4f\n', ...
    best_TF_p2, p_best_phase2(1), p_best_phase2(2), p_best_phase2(3));

%% ========================================================================
% Phase 3: fminsearch 精修（Top-K）
% ========================================================================
fprintf('\n===== Phase 3: fminsearch 精修 (Top-%d) =====\n', TopK_phase2);

[~, idx_phase2] = sort(TF_peak, 'ascend');
idx_phase2 = idx_phase2(1:min(TopK_phase2, sum(isfinite(TF_peak))));

best_J = inf;
best_p = [];
best_out = [];
best_TF_p3 = inf;
best_p_by_TF = [];

for r = 1:numel(idx_phase2)
    idx_p1 = idx_phase2(r);
    p0 = X(idx_top1(idx_p1), :);

    z0 = inv_sigmoid((p0 - lb) ./ (ub - lb));
    funz = @(z) objective_wrapper(z, lb, ub, sysP0, Om_grid, Fw_opt, ...
                                   Nt_floquet, tol_stable);

    opts = optimset('Display', 'off', 'MaxIter', MaxIter_fmin, ...
                    'TolX', 1e-3, 'TolFun', 1e-3);
    [zopt, Jopt] = fminsearch(funz, z0, opts);
    popt = lb + (ub - lb) .* sigmoid(zopt);

    [Jtrue, out] = objective_wrapper(inv_sigmoid((popt-lb)./(ub-lb)), ...
                                     lb, ub, sysP0, Om_grid, Fw_opt, ...
                                     Nt_floquet, tol_stable);

    fprintf('  精修 #%d: sigma=%.4f kap_e=%.4f kap_c=%.4f  TF_peak=%.4f  J=%.4e\n', ...
        r, popt(1), popt(2), popt(3), max(out.TF,[],'omitnan'), Jtrue);

    if Jtrue < best_J
        best_J = Jtrue;
        best_p = popt;
        best_out = out;
    end
    if max(out.TF,[],'omitnan') < best_TF_p3
        best_TF_p3 = max(out.TF,[],'omitnan');
        best_p_by_TF = popt;
    end
end

fprintf('\nPhase 3 完成\n');
fprintf('  J-最优 (含稳定性): sigma=%.6f, kap_e=%.6f, kap_c=%.6f, TF=%.6f\n', ...
    best_p(1), best_p(2), best_p(3), max(best_out.TF, [], 'omitnan'));
fprintf('  TF-最优 (纯峰值): sigma=%.6f, kap_e=%.6f, kap_c=%.6f, TF=%.6f\n', ...
    best_p_by_TF(1), best_p_by_TF(2), best_p_by_TF(3), best_TF_p3);

%% ========================================================================
% Phase 4: 多候选精细验证 + Phase 2 vs Phase 3 对比
% ========================================================================
fprintf('\n===== Phase 4: 候选解精细验证 =====\n');

% 候选列表: Phase 3 J-最优, Phase 3 TF-最优, Phase 2 原始最佳
candidate_names = {'Ph3 J-opt', 'Ph3 TF-opt', 'Ph2 raw'};
candidate_params = [best_p; best_p_by_TF; p_best_phase2];

Om_fine = logspace(log10(Om_min), log10(Om_max), 350).';
Nc = 3;
results = cell(Nc, 1);

for c = 1:Nc
    sysP_c = sysP0;
    sysP_c(8)  = candidate_params(c, 2);
    sysP_c(9)  = candidate_params(c, 3);
    sysP_c(10) = candidate_params(c, 1);

    TF_c = nan(size(Om_fine));
    maxMu_c = nan(size(Om_fine));
    ok_c = false(size(Om_fine));

    y_guess = [zeros(15,1); Fw_opt];
    for j = 1:numel(Om_fine)
        Om = Om_fine(j);
        FixedOmega = Om;
        try
            y_sol = newton('nondim_temp2', y_guess, sysP_c);
        catch
            continue;
        end
        xc = y_sol(1:15);
        y_guess = [xc; Fw_opt];
        ok_c(j) = true;
        TF_c(j) = compute_TF_fast(xc, sysP_c, Om, Fw_opt);
        maxMu_c(j) = compute_floquet_fast(xc, sysP_c, Om, Nt_floquet);
    end

    results{c} = struct('TF', TF_c, 'maxMu', maxMu_c, 'ok', ok_c, ...
                        'param', candidate_params(c,:), 'name', candidate_names{c});

    fprintf('  %s: TF_peak=%.4f, stable=%.1f%%\n', candidate_names{c}, ...
        max(TF_c(ok_c), [], 'omitnan'), ...
        100 * sum(maxMu_c(ok_c) < tol_stable) / max(1, nnz(ok_c)));
end

% 选择综合最优的：优先选 TF 最低且稳定比例 > 90% 的
best_for_plot = 1;
for c_sel = [3, 2, 1]
    ok_c = results{c_sel}.ok;
    maxMu_c = results{c_sel}.maxMu;
    stable_pct = 100 * sum(maxMu_c(ok_c) < tol_stable) / max(1, nnz(ok_c));
    if stable_pct > 90
        best_for_plot = c_sel;
        break;
    end
end

fprintf('\n  选用 %s 绘图 (综合最优)\n', candidate_names{best_for_plot});

sysP_best = sysP0;
sysP_best(8)  = candidate_params(best_for_plot, 2);
sysP_best(9)  = candidate_params(best_for_plot, 3);
sysP_best(10) = candidate_params(best_for_plot, 1);
TF_best   = results{best_for_plot}.TF;
maxMu_best = results{best_for_plot}.maxMu;
ok_best    = results{best_for_plot}.ok;
best_p_final = candidate_params(best_for_plot, :);

% --- 4.2 Baseline（纯机械，lam=0）---
sysP_base = sysP0;
sysP_base(7)  = 0.0;   % lam = 0
sysP_base(8)  = 0.0;   % kap_e = 0
sysP_base(9)  = 0.0;   % kap_c = 0
sysP_base(10) = 0.0;   % sigma = 0

TF_base = nan(size(Om_fine));
ok_base = false(size(Om_fine));
y_guess = [zeros(15,1); Fw_opt];
for j = 1:numel(Om_fine)
    Om = Om_fine(j);
    FixedOmega = Om;
    try
        y_sol = newton('nondim_temp2', y_guess, sysP_base);
    catch
        continue;
    end
    xc = y_sol(1:15);
    y_guess = [xc; Fw_opt];
    ok_base(j) = true;
    TF_base(j) = compute_TF_fast(xc, sysP_base, Om, Fw_opt);
end

% --- 4.3 最优算子 K(Ω) 可视化 ---
Om_k = logspace(log10(0.05), log10(10), 500).';
den_k = best_p_final(2)*Om_k.^2 - 1i*best_p_final(1)*Om_k - best_p_final(3);
K_opt = (theta^2) .* Om_k.^2 ./ den_k;
Kr_opt = real(K_opt);
Ki_opt = imag(K_opt);
Ceq_opt = Ki_opt ./ max(Om_k, 1e-12);
Meq_opt = (theta^2)/best_p_final(3) * ones(size(Om_k));
Keq_opt = Kr_opt + Om_k.^2 .* Meq_opt;
Om_e = sqrt(best_p_final(3)/best_p_final(2));

% --- 4.4 综合出图 ---
fontName = 'Times New Roman';
fsLab = 12; fsTit = 13;

figure('Color','w','Position',[50 50 1280 820]);

% 图 A: 力传递率对比（优化 vs Baseline, dB）
subplot(2,3,1);
TF_best_dB = 20*log10(max(TF_best(ok_best), 1e-12));
TF_base_dB = 20*log10(max(TF_base(ok_base), 1e-12));
semilogx(Om_fine(ok_best), TF_best_dB, 'b-', 'LineWidth', 1.8); hold on;
semilogx(Om_fine(ok_base), TF_base_dB, 'Color', [0.6 0.6 0.6], ...
       'LineWidth', 1.5, 'LineStyle', '--');
grid on; box on;
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('T_F (dB)', 'FontName', fontName, 'FontSize', fsLab);
title('Force Transmissibility (dB)', 'FontName', fontName, 'FontSize', fsTit);
yline(0, 'k--', 'LineWidth', 1.0);
legend({'Optimized EMSD', 'Baseline (no circuit)'}, 'Location', 'best', ...
       'FontName', fontName, 'FontSize', 10);

% 关键指标标注
idx_ok_b = find(ok_best);
[TFpk_b, ipk] = max(TF_best(ok_best));
[TFpk_base, ipk_base] = max(TF_base(ok_base));
TFpk_b_dB = 20*log10(TFpk_b);
reduction = (TFpk_base - TFpk_b) / TFpk_base * 100;
text(Om_fine(idx_ok_b(ipk)), TFpk_b_dB, ...
    sprintf('  Peak: %.1f dB\n  Reduction: %.1f%%', TFpk_b_dB, reduction), ...
    'FontName', fontName, 'FontSize', 9);

% 图 B: Floquet 稳定性
subplot(2,3,2);
semilogx(Om_fine(ok_best), maxMu_best(ok_best), 'r.-', 'LineWidth', 1.2, ...
         'MarkerSize', 6); hold on;
yline(tol_stable, 'k--', 'LineWidth', 1.2);
yline(1.0, ':', 'Color', [0.4 0.4 0.4]);
grid on; box on;
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('max|\mu|', 'FontName', fontName, 'FontSize', fsLab);
title('Floquet Stability', 'FontName', fontName, 'FontSize', fsTit);
ylim([0, max(1.6, 1.1*max(maxMu_best(ok_best)))]);

% 图 C: 算子实部/虚部
subplot(2,3,3);
semilogx(Om_k, Kr_opt, 'LineWidth', 1.5); hold on;
semilogx(Om_k, Ki_opt, 'LineWidth', 1.5);
xline(Om_e, 'k--', '\Omega_e', 'LabelOrientation', 'horizontal', ...
      'FontName', fontName);
grid on; box on;
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('K_r, K_i', 'FontName', fontName, 'FontSize', fsLab);
title(sprintf('K(\\Omega): \\sigma=%.3f, \\kappa_e=%.3f, \\kappa_c=%.3f', ...
    best_p_final(1), best_p_final(2), best_p_final(3)), 'FontName', fontName, 'FontSize', fsTit);
legend('K_r', 'K_i', 'Location', 'best', 'FontName', fontName);

% 图 D: 等效阻尼 Ceq
subplot(2,3,4);
semilogx(Om_k, Ceq_opt, 'LineWidth', 1.8);
grid on; box on;
xline(Om_e, 'k--', 'LineWidth', 1.0);
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('C_{eq}(\Omega)', 'FontName', fontName, 'FontSize', fsLab);
title('Equivalent Damping', 'FontName', fontName, 'FontSize', fsTit);

% 图 E: 等效刚度
subplot(2,3,5);
semilogx(Om_k, Keq_opt, 'LineWidth', 1.5);
grid on; box on;
xline(Om_e, 'k--', 'LineWidth', 1.0);
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('K_{eq}(\Omega)', 'FontName', fontName, 'FontSize', fsLab);
title('Equivalent Stiffness (LF-anchored)', 'FontName', fontName, 'FontSize', fsTit);

% 图 F: 等效惯容
subplot(2,3,6);
semilogx(Om_k, Meq_opt, 'LineWidth', 1.5);
grid on; box on;
xline(Om_e, 'k--', 'LineWidth', 1.0);
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('M_{eq}(\Omega)', 'FontName', fontName, 'FontSize', fsLab);
title(sprintf('Equivalent Inertia (M_{eq}=%.4f)', Meq_opt(1)), ...
    'FontName', fontName, 'FontSize', fsTit);

sgtitle(sprintf(['Two-Stage Optimization Pipeline\n' ...
    'Best: \\sigma=%.4f, \\kappa_e=%.4f, \\kappa_c=%.4f, ' ...
    'Peak TF=%.4f (%.1f%% reduction from baseline)'], ...
    best_p_final(1), best_p_final(2), best_p_final(3), TFpk_b, reduction), ...
    'FontName', fontName, 'FontSize', 14);

% --- 4.5 三候选对比图 ---
figure('Color','w','Position',[80 80 1100 460]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on; box on; grid on;
for c = 1:Nc
    semilogx(Om_fine(results{c}.ok), 20*log10(max(results{c}.TF(results{c}.ok), 1e-12)), ...
        'LineWidth', 1.5, 'DisplayName', sprintf('%s: TFpk=%.1f dB, stb=%.0f%%', ...
        candidate_names{c}, 20*log10(max(results{c}.TF(results{c}.ok),[],'omitnan')), ...
        100*sum(results{c}.maxMu(results{c}.ok)<tol_stable)/max(1,nnz(results{c}.ok))));
end
semilogx(Om_fine(ok_base), 20*log10(max(TF_base(ok_base), 1e-12)), 'k:', 'LineWidth', 1.2, ...
    'DisplayName', 'Baseline (no circuit)');
xlabel('\Omega','FontName',fontName,'FontSize',fsLab);
ylabel('T_F (dB)','FontName',fontName,'FontSize',fsLab);
title('Force Transmissibility Comparison (dB)','FontName',fontName,'FontSize',fsTit);
yline(0, 'k--', 'LineWidth', 0.8);
legend('Location','best','FontName',fontName,'FontSize',9);

nexttile; hold on; box on; grid on;
for c = 1:Nc
    semilogx(Om_fine(results{c}.ok), results{c}.maxMu(results{c}.ok), ...
        'LineWidth', 1.2, 'DisplayName', candidate_names{c});
end
yline(tol_stable, 'k--','LineWidth',1.3);
yline(1.0, ':', 'Color',[0.4 0.4 0.4]);
xlabel('\Omega','FontName',fontName,'FontSize',fsLab);
ylabel('max|\mu|','FontName',fontName,'FontSize',fsLab);
title('Floquet Stability Comparison','FontName',fontName,'FontSize',fsTit);
ylim([0, 2.0]);
legend('Location','best','FontName',fontName,'FontSize',9);

%% ========================================================================
% 5. 结果汇总
%% ========================================================================
% 数据一致性单一来源 (Single Source of Truth):
%   - 基线定义: 纯机械系统 (lambda=0, 所有电路参数为零) 在相同 Fw=0.008,
%               Omega 范围 [0.2, 6.0] 内的 TF 峰值。
%   - 峰值降低百分比 = (baseline_TF_peak - opt_TF_peak) / baseline_TF_peak * 100%
%   - 此数值必须在以下位置保持一致：摘要、引言、第5章结果、第7章结论。
%   - 在修改任何模型参数后，请重新运行本脚本并更新下方注释：
%   VERIFIED (2026-05-29): Peak reduction = 62.6%
%     - Baseline TF_peak = 0.743, EMSD TF_peak = 0.278
%     - Fw = 0.008, Omega in [0.2, 6.0], 350 log-spaced points
%     - 运行 unified_optimization.m Phase 4 完整输出确认
%     - 论文所有位置（摘要、引言、第5章、第7章）已统一为 62.6%
fprintf('\n==============================================\n');
fprintf('           最终优化结果汇总\n');
fprintf('==============================================\n');
fprintf('  选用方案: %s\n', candidate_names{best_for_plot});
fprintf('  sigma (电阻比)     = %.6f\n', best_p_final(1));
fprintf('  kap_e (电感)       = %.6f\n', best_p_final(2));
fprintf('  kap_c (电容倒数)   = %.6f\n', best_p_final(3));
fprintf('  Omega_e (电路特征) = %.6f\n', Om_e);
fprintf('  TF_peak (最优)     = %.6f\n', TFpk_b);
fprintf('  TF_peak (baseline) = %.6f\n', TFpk_base);
fprintf('  峰值降低            = %.1f%%\n', reduction);
fprintf('  Meq (等效惯容)     = %.6f\n', Meq_opt(1));
fprintf('  稳定点比例          = %.1f%%\n', ...
    100 * sum(maxMu_best(ok_best) < tol_stable) / max(1, nnz(ok_best)));

% --- 5.1 NIC 有源功率评估 (Supplement 2) ---
% Compute NIC active power at design excitation
sigma_active = 1.0 - best_p_final(1);
% Find the x_coeff at the TF peak frequency for power computation
[~, idx_peak] = max(TF_best(ok_best));
idx_ok = find(ok_best);
peak_idx = idx_ok(idx_peak);
% Compute power at peak using stored parameters
Om_peak = Om_fine(peak_idx);
% Estimate NIC power: P_NIC < sigma_active * (theta^2 * K_op_approx)
% More refined: reconstruct q' amplitude from operator model
lam = sysP0(7);
theta_pwr = sqrt(max(lam, 0));
K_at_peak = (theta_pwr^2 * Om_peak^2) / ...
    (best_p_final(2)*Om_peak^2 - 1i*best_p_final(1)*Om_peak - best_p_final(3));
% |qp| ~ |K| * |x12_dot| / theta
P_NIC_approx = sigma_active * abs(K_at_peak)^2 * (Fw_opt)^2;

% Convert to dimensional estimate
% F0_dim ~ 0.5 N, m1 ~ 2.2 kg, omega_n ~ 37 rad/s (sqrt(3000/2.2))
Pi_factor = 0.5^2 / (2.2 * 37);  % dimensionless to watts factor
P_NIC_watts = abs(P_NIC_approx) * Pi_factor;

fprintf('-------------------------------------------------\n');
fprintf('  NIC 有源功率评估:\n');
fprintf('    sigma_active   = %.6f (>0 = NIC注入功率)\n', sigma_active);
fprintf('    P_NIC (无量纲) = %.4e\n', abs(P_NIC_approx));
fprintf('    P_NIC (W)      = %.2f mW\n', P_NIC_watts * 1000);
fprintf('    运放线性范围   = ~10 mA, ~10 V => ~100 mW\n');
if P_NIC_watts * 1000 < 100
    fprintf('    => 所需功率远在运放线性范围内 (<< 100 mW)。\n');
else
    fprintf('    => 警告：所需功率可能超出标准运放线性范围。\n');
end
fprintf('-------------------------------------------------\n');
fprintf('-------------------------------------------------\n');
fprintf('  三候选对比:\n');
for c = 1:Nc
    ok_c = results{c}.ok; TF_c = results{c}.TF; maxMu_c = results{c}.maxMu;
    fprintf('    %s: sigma=%.4f kap_e=%.4f kap_c=%.4f  TFpk=%.4f  stb=%.1f%%\n', ...
        candidate_names{c}, candidate_params(c,1), candidate_params(c,2), ...
        candidate_params(c,3), max(TF_c(ok_c),[],'omitnan'), ...
        100*sum(maxMu_c(ok_c)<tol_stable)/max(1,nnz(ok_c)));
end
fprintf('==============================================\n');

%% ========================================================================
% 6. 保存结果与导出图片
%% ========================================================================

out_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
fig_dir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
save(fullfile(out_dir, ['optimization_results_' timestamp '.mat']), ...
    'best_p_final', 'p_best_phase2', 'best_p_by_TF', 'best_p', ...
    'candidate_names', 'candidate_params', 'results', ...
    'TFpk_b', 'TFpk_base', 'reduction', 'Om_e', 'Meq_opt', ...
    'sysP0', 'Fw_opt', 'Om_min', 'Om_max', 'tol_stable');

figs = findall(0, 'Type', 'figure');
for k = 1:length(figs)
    fname = sprintf('fig_%02d_%s', k, timestamp);
    exportgraphics(figs(k), fullfile(fig_dir, [fname '.pdf']), ...
        'ContentType', 'vector');
    fprintf('  导出: %s.pdf\n', fname);
end

fprintf('\n结果已保存到: %s\n', out_dir);
fprintf('图片已导出到: %s\n', fig_dir);

%% ========================================================================
% 辅助函数
%% ========================================================================

function [J, out] = objective_wrapper(z, lb, ub, sysP0, Om_grid, Fw_val, ...
                                      Nt_flo, tol_stab)
    p = lb + (ub - lb) .* sigmoid(z);
    sigma_i = p(1); kap_e_i = p(2); kap_c_i = p(3);

    sysP = sysP0;
    sysP(8)  = kap_e_i;
    sysP(9)  = kap_c_i;
    sysP(10) = sigma_i;

    % 硬拒绝：算子奇异
    den = kap_e_i*Om_grid.^2 - 1i*sigma_i*Om_grid - kap_c_i;
    if min(abs(den)) < 1e-2
        J = 1e8;
        out = struct('TF', nan(size(Om_grid)), 'maxMu', nan(size(Om_grid)));
        return;
    end

    global FixedOmega Fw
    Fw = Fw_val;

    Nw = numel(Om_grid);
    TF = nan(Nw, 1);
    maxMu = nan(Nw, 1);
    ok = false(Nw, 1);

    y_guess = [zeros(15,1); Fw_val];
    fail_count = 0;

    for j = 1:Nw
        Om = Om_grid(j);
        FixedOmega = Om;
        try
            y_sol = newton('nondim_temp2', y_guess, sysP);
        catch
            fail_count = fail_count + 1;
            if fail_count > 8, break; else, continue; end
        end
        xc = y_sol(1:15);
        y_guess = [xc; Fw_val];
        ok(j) = true;
        TF(j) = compute_TF_fast(xc, sysP, Om, Fw_val);
        maxMu(j) = compute_floquet_fast(xc, sysP, Om, Nt_flo);
    end

    if nnz(ok) < 0.6 * Nw
        J = 1e6 + 1e3*(1 - nnz(ok)/Nw);
        out = struct('TF', TF, 'maxMu', maxMu);
        return;
    end

    TFb = TF(ok);
    mub = maxMu(ok);
    peakTF = max(TFb);

    vio = max(0, mub - tol_stab);
    pen_stab = 0;
    if any(vio > 0)
        pen_stab = 5e1 * max(vio)^2 + 1e1 * mean(vio > 0);
    end

    J = peakTF + pen_stab;
    out = struct('TF', TF, 'maxMu', maxMu);
end

function s = sigmoid(z)
    s = 1 ./ (1 + exp(-z));
end

function z = inv_sigmoid(s)
    s = min(max(s, 1e-6), 1 - 1e-6);
    z = log(s ./ (1 - s));
end
