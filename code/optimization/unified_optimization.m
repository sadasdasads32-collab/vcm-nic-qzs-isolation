%% unified_optimization.m
% =========================================================================
% 两级优化流水线：HBM快筛 → 系统精修 → 稳定性验证
%
% Phase 1: 系统层 HBM/弧长延拓快筛（800候选，Budget=1500）
%          → 直接 HBM 评估 TF 峰值 + Floquet 抽查，选 Top-30
% Phase 2: 系统层 FRF 验证（HBM+弧长延拓，Top-30 候选）
%          → 高精度 FRF (Budget=3000)，精选 Top-8
% Phase 3: 系统层 fminsearch 精修（Top-8 候选局部优化 + 提前终止）
%          → 在真实 HBM 模型上微调参数，收敛后精细验证
% Phase 4: Floquet 稳定性全扫描 + Wang BG 模型基线 + 算子解释 + 出版级对比图
%
% v3.0 改进:
%   - NIC 负参数搜索: sigma ∈ [-3, 3], kap_c ∈ [-3, 3]
%   - Phase 1 改为直接 HBM 快筛 (跳过算子预筛，避免负参数评分失效)
%   - 新增 Wang 2017 BG 模型基线对比
%   - 新增算子 K(Ω) 物理解释 (6-panel 分解图)
%   - 新增 Wang 风格参数扫掠图
%   - 新增 4 项 Wang 性能指标对比表
% =========================================================================

clc; clear; close all;
init_path();
% 抑制弧长延拓中近奇异矩阵的冗余警告（不影响结果，仅减少 I/O）
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');

%% ========================================================================
% 0. 用户可调参数
% ========================================================================

% --- 机械参数（与 Wang(2017) BG 模型一致）---
mu   = 0.2;         % 质量比
beta = 2.0;         % 下层线性刚度比
K1   = 1.0;         % 上层水平弹簧比
K2   = 0.0;         % 下层水平弹簧比
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
Fw_opt  = 0.005;    % 优化用激励幅值
Om_min  = 0.2;
Om_max  = 6.0;

% --- Phase 1: HBM 快筛 ---
Nsamp_phase1  = 800;      % LHS 采样数（直接 HBM 评估）
TopN_phase1   = 30;       % 进入 Phase 2 的候选数
Budget_p1     = 1500;     % Phase 1 弧长预算（轻量化）

% --- Phase 2: FRF 验证 ---
TopK_phase2   = 8;        % 进入 Phase 3 精修的候选数
Budget_p2     = 3000;     % Phase 2 弧长预算

% --- Phase 3: fminsearch ---
MaxIter_fmin  = 80;       % fminsearch 最大迭代

% --- 搜索边界 [sigma, kap_e, kap_c] ---
% NIC 允许 sigma (负阻) 和 kap_c (负容抗) 取负值
lb = [-3.00, 0.02, -3.00];
ub = [ 3.00, 3.00,  3.00];

% --- Floquet 设置 ---
Nt_floquet  = 600;
tol_stable  = 1.002;

% --- 并行计算设置 ---
% 注意: 部分MATLAB版本不支持parfor内声明global变量
% 若遇到"在parfor循环中不支持全局变量"错误，请设置为false
use_parallel = false;
try
    pool = gcp('nocreate');
    if isempty(pool) && use_parallel
        pool = parpool('local');
    end
    if ~isempty(pool)
        fprintf('  并行池就绪, 工作进程数: %d\n', pool.NumWorkers);
    end
catch ME
    fprintf('  并行池不可用 (%s), 回退到串行模式\n', ME.message);
    use_parallel = false;
end

rng(42);
fprintf('========== 优化流水线 v3.0 启动 ==========\n');
fprintf('lam=%.4f (theta=%.4f), Fw=%.4f\n', lam_phys, theta, Fw_opt);
fprintf('搜索边界: sigma∈[%.2f,%.2f], kap_e∈[%.2f,%.2f], kap_c∈[%.2f,%.2f]\n\n', ...
    lb(1), ub(1), lb(2), ub(2), lb(3), ub(3));

%% ========================================================================
% Phase 1: 系统层 HBM/弧长延拓快筛（800 候选，Budget=1500）
% ========================================================================
fprintf('===== Phase 1: HBM/弧长延拓快筛 (%d 样本, Budget=%d) =====\n', ...
    Nsamp_phase1, Budget_p1);

X = lhsdesign(Nsamp_phase1, 3);
X = lb + X .* (ub - lb);

TF_peak_p1  = inf(Nsamp_phase1, 1);
Mu_peak_p1  = inf(Nsamp_phase1, 1);

global FixedOmega Fw ParamMin ParamMax
Fw = Fw_opt;
FixedOmega = [];
ParamMin = 0.05;
ParamMax = 10.5;

sysP0 = [be1, be2, mu, al1, gamma1, ze1, lam_phys, ...
         1.0, 0.2, 1.0, gamma2];

tic_phase1 = tic;
for i = 1:Nsamp_phase1
    sigma_i = X(i,1); kap_e_i = X(i,2); kap_c_i = X(i,3);

    % 跳过 kap_c 过零附近的奇异点 (1/kap_c → ∞)
    if abs(kap_c_i) < 1e-4
        continue;
    end

    sysP = sysP0;
    sysP(8)  = kap_e_i;
    sysP(9)  = kap_c_i;
    sysP(10) = sigma_i;

    % Phase 1 使用 arc_length_frf_robust 直接评估（无 fallback，避免冗长输出）
    try
        [Om, TF_dB, x_res] = arc_length_frf_robust(sysP, 10.0, ...
            'Fw', Fw_opt, 'Budget', Budget_p1, 'Verbose', false);
    catch
        continue;
    end
    if isempty(Om) || length(Om) < 10
        continue;
    end

    % 线性 TF 峰值
    TF_lin = 10.^(TF_dB/20);
    [peak_val, idx_peak] = max(TF_lin);
    TF_peak_p1(i) = peak_val;

    % Floquet 稳定性抽查 (低频/峰值/高频 3 点)
    Om_raw = x_res(16, :)';
    valid_x = Om_raw > 0 & isfinite(Om_raw);
    if nnz(valid_x) < 3
        continue;
    end
    x_res_f = x_res(:, valid_x);
    Om_peak_val = Om(idx_peak);
    [~, idx_peak_x] = min(abs(x_res_f(16, :)' - Om_peak_val));
    check_idx = unique([1, idx_peak_x, size(x_res_f, 2)]);
    mu_vals = zeros(1, length(check_idx));
    for k = 1:length(check_idx)
        xc_k = x_res_f(1:15, check_idx(k));
        Om_k = x_res_f(16, check_idx(k));
        mu_vals(k) = compute_floquet_fast(xc_k, sysP, Om_k, Nt_floquet);
    end
    Mu_peak_p1(i) = max(mu_vals);

    if mod(i, 80) == 0
        fprintf('  Phase 1: %3d/%d done, best TF=%.4f\n', ...
            i, Nsamp_phase1, min(TF_peak_p1(1:i)));
    end
end
t_phase1 = toc(tic_phase1);

% 按 TF 峰值升序排列 (越小越好 == 更好的隔离)
[~, idx_sort] = sort(TF_peak_p1, 'ascend');
idx_valid = idx_sort(isfinite(TF_peak_p1(idx_sort)));
TopN_actual = min(TopN_phase1, numel(idx_valid));
idx_top1 = idx_valid(1:TopN_actual);

fprintf('Phase 1 完成 (%.1f s), Top-%d 候选已选出\n', t_phase1, TopN_actual);
fprintf('  最佳 TF_peak = %.4f\n', TF_peak_p1(idx_top1(1)));
fprintf('  参数: sigma=%.4f, kap_e=%.4f, kap_c=%.4f\n', ...
    X(idx_top1(1),1), X(idx_top1(1),2), X(idx_top1(1),3));

% 统计负参数候选在 Top-30 中的比例
neg_sigma_count = sum(X(idx_top1, 1) < 0);
neg_kapc_count  = sum(X(idx_top1, 3) < 0);
fprintf('  Top-%d 中 sigma<0: %d 个, kap_c<0: %d 个\n\n', ...
    TopN_actual, neg_sigma_count, neg_kapc_count);

%% ========================================================================
% Phase 2: 系统层 FRF 验证（Top-N，弧长延拓，可选并行）
% ========================================================================
fprintf('===== Phase 2: 系统层 FRF 验证 (Top-%d, Budget=%d', TopN_actual, Budget_p2);
if use_parallel, fprintf(', 并行'); end
fprintf(') =====\n');

cand_X = X(idx_top1(1:TopN_actual), :);
TF_peak  = inf(TopN_actual, 1);
Mu_peak  = inf(TopN_actual, 1);

tic_phase2 = tic;

if use_parallel
    % ==== 并行路径 ====
    Fw_p2  = Fw_opt;
    Nt_p2  = Nt_floquet;
    parfor r = 1:TopN_actual
        init_path();  % 确保 worker 上有库路径
        global Fw ParamMin ParamMax FixedOmega
        Fw = Fw_p2;
        FixedOmega = [];
        ParamMin = 0.05;
        ParamMax = 10.5;

        sysP_r = sysP0;
        sysP_r(8)  = cand_X(r, 2);
        sysP_r(9)  = cand_X(r, 3);
        sysP_r(10) = cand_X(r, 1);

        TF_pk = 1e6;
        Mu_pk = 1e6;

        try
            [Om, TF_dB, x_res] = call_robust_frf(sysP_r, Fw_p2, Budget_p2);
            if ~isempty(Om) && ~isempty(TF_dB)
                TF_lin = 10.^(TF_dB/20);
                [TF_pk, idx_peak] = max(TF_lin);

                % Floquet 稳定性抽查 (低频/峰值/高频)
                Om_raw = x_res(16, :)';
                valid_x = Om_raw > 0 & isfinite(Om_raw);
                x_res_f = x_res(:, valid_x);
                Om_peak_val = Om(idx_peak);
                [~, idx_peak_x] = min(abs(x_res_f(16, :)' - Om_peak_val));
                check_idx = unique([1, idx_peak_x, size(x_res_f, 2)]);
                mu_vals = zeros(1, length(check_idx));
                for k = 1:length(check_idx)
                    xc_k = x_res_f(1:15, check_idx(k));
                    Om_k = x_res_f(16, check_idx(k));
                    mu_vals(k) = compute_floquet_fast(xc_k, sysP_r, Om_k, Nt_p2);
                end
                Mu_pk = max(mu_vals);
            end
        catch
        end
        TF_peak(r) = TF_pk;
        Mu_peak(r) = Mu_pk;
    end
else
    % ==== 串行路径 (含进度输出) ====
    for r = 1:TopN_actual
        sigma_i = cand_X(r, 1); kap_e_i = cand_X(r, 2); kap_c_i = cand_X(r, 3);

        sysP = sysP0;
        sysP(8)  = kap_e_i;
        sysP(9)  = kap_c_i;
        sysP(10) = sigma_i;

        % --- 弧长延拓 FRF ---
        try
            [Om, TF_dB, x_res] = call_robust_frf(sysP, Fw_opt, Budget_p2);
        catch
            TF_peak(r) = 1e6;
            Mu_peak(r) = 1e6;
            if mod(r, 5) == 0
                fprintf('  Phase 2: %2d/%d done (arc-length failed)\n', r, TopN_actual);
            end
            continue;
        end

        if isempty(Om) || isempty(TF_dB)
            TF_peak(r) = 1e6;
            Mu_peak(r) = 1e6;
            continue;
        end

        % 线性 TF 峰值
        TF_lin = 10.^(TF_dB/20);
        [peak_TF_val, idx_peak] = max(TF_lin);
        TF_peak(r) = peak_TF_val;

        % --- Floquet 稳定性抽查（低频/峰值/高频 3 个去重点）---
        Om_raw = x_res(16, :)';
        valid_x = Om_raw > 0 & isfinite(Om_raw);
        x_res_f = x_res(:, valid_x);
        Om_peak_val = Om(idx_peak);
        [~, idx_peak_x] = min(abs(x_res_f(16, :)' - Om_peak_val));

        check_idx = unique([1, idx_peak_x, size(x_res_f, 2)]);
        mu_vals = zeros(1, length(check_idx));
        for k = 1:length(check_idx)
            xc_k = x_res_f(1:15, check_idx(k));
            Om_k = x_res_f(16, check_idx(k));
            mu_vals(k) = compute_floquet_fast(xc_k, sysP, Om_k, Nt_floquet);
        end
        Mu_peak(r) = max(mu_vals);

        if mod(r, 5) == 0
            fprintf('  Phase 2: %2d/%d done, best TF_peak=%.3f\n', ...
                r, TopN_actual, min(TF_peak(1:r)));
        end
    end
end
t_phase2 = toc(tic_phase2);

% 保存 Phase 2 原始最佳（用于后续对比）
[best_TF_p2, idx_best_p2] = min(TF_peak);
p_best_phase2 = X(idx_top1(idx_best_p2), :);

fprintf('Phase 2 完成 (%.1f s)\n', t_phase2);
fprintf('  最佳 TF_peak = %.4f @ sigma=%.4f kap_e=%.4f kap_c=%.4f\n\n', ...
    best_TF_p2, p_best_phase2(1), p_best_phase2(2), p_best_phase2(3));

%% ========================================================================
% Phase 3: fminsearch 精修（Top-K，提前终止 + 精细验证）
% ========================================================================
fprintf('===== Phase 3: fminsearch 精修 (Top-%d) =====\n', TopK_phase2);

[~, idx_phase2] = sort(TF_peak, 'ascend');
idx_phase2 = idx_phase2(1:min(TopK_phase2, sum(isfinite(TF_peak))));

best_J = inf;
best_p = [];
best_out = [];
best_TF_p3 = inf;
best_p_by_TF = [];
MAX_STALL = 15;  % 提前终止阈值: 连续 N 次无 1% 改善

for r = 1:numel(idx_phase2)
    idx_p1 = idx_phase2(r);
    p0 = X(idx_top1(idx_p1), :);

    z0 = inv_sigmoid((p0 - lb) ./ (ub - lb));
    funz = @(z) objective_wrapper(z, lb, ub, sysP0, Fw_opt, ...
                                   Nt_floquet, tol_stable);

    % fminsearch 带提前终止
    early_stop = make_early_stop(MAX_STALL);
    opts = optimset('Display', 'off', 'MaxIter', MaxIter_fmin, ...
                    'TolX', 1e-3, 'TolFun', 1e-3, ...
                    'OutputFcn', early_stop);
    [zopt, Jopt] = fminsearch(funz, z0, opts);
    popt = lb + (ub - lb) .* sigmoid(zopt);

    % 真实评估 (使用 objective_wrapper 的精细弧长)
    [Jtrue, out] = objective_wrapper(inv_sigmoid((popt-lb)./(ub-lb)), ...
                                     lb, ub, sysP0, Fw_opt, ...
                                     Nt_floquet, tol_stable);

    % --- 收敛后精细验证: 3000步弧长 vs 目标函数弧长 ---
    TF_fine_peak = out.TF_peak;
    try
        sysP_val = sysP0;
        sysP_val(8)  = popt(2);
        sysP_val(9)  = popt(3);
        sysP_val(10) = popt(1);
        global ParamMin ParamMax
        ParamMin = 0.05; ParamMax = 10.5;
        [~, TF_dB_fine] = call_robust_frf(sysP_val, Fw_opt, 3000);
        if ~isempty(TF_dB_fine)
            TF_fine_peak = max(10.^(TF_dB_fine/20));
        end
    catch
    end

    disc = abs(TF_fine_peak - out.TF_peak) / max(TF_fine_peak, 1e-6);

    if disc > 0.20
        fprintf('  #%d: sig=%.4f kap_e=%.4f kap_c=%.4f  coarse_TF=%.4f fine_TF=%.4f ***Delta=%.0f%%***\n', ...
            r, popt(1), popt(2), popt(3), out.TF_peak, TF_fine_peak, disc*100);
        pen_stab = Jtrue - out.TF_peak;  % 隔离稳定性惩罚 (在覆盖 out.TF_peak 之前)
        out.TF_peak = TF_fine_peak;
        Jtrue = TF_fine_peak + pen_stab;
    else
        fprintf('  #%d: sig=%.4f kap_e=%.4f kap_c=%.4f  TF=%.4f  J=%.4e\n', ...
            r, popt(1), popt(2), popt(3), out.TF_peak, Jtrue);
    end

    if Jtrue < best_J
        best_J = Jtrue;
        best_p = popt;
        best_out = out;
    end
    if out.TF_peak < best_TF_p3
        best_TF_p3 = out.TF_peak;
        best_p_by_TF = popt;
    end
end

fprintf('\nPhase 3 完成\n');
fprintf('  J-最优 (含稳定性): sigma=%.6f, kap_e=%.6f, kap_c=%.6f, TF=%.6f\n', ...
    best_p(1), best_p(2), best_p(3), best_out.TF_peak);
fprintf('  TF-最优 (纯峰值): sigma=%.6f, kap_e=%.6f, kap_c=%.6f, TF=%.6f\n', ...
    best_p_by_TF(1), best_p_by_TF(2), best_p_by_TF(3), best_TF_p3);

%% ========================================================================
% Phase 4: 多候选精细验证 + Wang BG 基线 + 算子解释 + 出版级图表
% ========================================================================
fprintf('\n===== Phase 4: 候选解精细验证 + Wang BG 基线 =====\n');

% 候选列表: Phase 3 J-最优, Phase 3 TF-最优, Phase 2 原始最佳
candidate_names = {'Ph3 J-opt', 'Ph3 TF-opt', 'Ph2 raw'};
candidate_params = [best_p; best_p_by_TF; p_best_phase2];

% Phase 4 使用 robust 弧长延拓（上扫 + 下扫回退）
global ParamMin ParamMax
ParamMin = 0.05;
ParamMax = 10.5;

Nc = 3;
results = cell(Nc, 1);
x_res_all = cell(Nc, 1);  % 保存完整 x_res 用于 Wang 指标计算

for c = 1:Nc
    sysP_c = sysP0;
    sysP_c(8)  = candidate_params(c, 2);
    sysP_c(9)  = candidate_params(c, 3);
    sysP_c(10) = candidate_params(c, 1);

    % --- 弧长延拓 FRF (robust: 上扫优先, 下扫回退) ---
    try
        [Om_c, TF_dB_c, x_res_c] = call_robust_frf(sysP_c, Fw_opt, 5000);
    catch
        results{c} = struct('TF', nan, 'maxMu', nan, 'ok', false, ...
            'Om', nan, 'Om_floq', nan, ...
            'param', candidate_params(c,:), 'name', candidate_names{c});
        x_res_all{c} = [];
        fprintf('  %s: arc-length failed\n', candidate_names{c});
        continue;
    end

    if isempty(Om_c)
        results{c} = struct('TF', nan, 'maxMu', nan, 'ok', false, ...
            'Om', nan, 'Om_floq', nan, ...
            'param', candidate_params(c,:), 'name', candidate_names{c});
        x_res_all{c} = [];
        continue;
    end

    x_res_all{c} = x_res_c;

    % 线性 TF
    TF_c_lin = 10.^(TF_dB_c/20);
    ok_c = true(size(Om_c));

    % --- Floquet 在 ~200 个自适应分布频点采样 ---
    N_floq_pts = min(200, length(Om_c));
    floq_idx = unique(round(linspace(1, length(Om_c), N_floq_pts)));
    N_pts = length(floq_idx);
    maxMu_c = nan(N_pts, 1);
    Om_floq = Om_c(floq_idx);

    for k = 1:N_pts
        xc_k = x_res_c(1:15, floq_idx(k));
        Om_k = x_res_c(16, floq_idx(k));
        maxMu_c(k) = compute_floquet_fast(xc_k, sysP_c, Om_k, Nt_floquet);
    end

    results{c} = struct('TF', TF_c_lin, 'TF_dB', TF_dB_c, ...
        'maxMu', maxMu_c, 'ok', ok_c, ...
        'Om', Om_c, 'Om_floq', Om_floq, ...
        'x_res', x_res_c, ...
        'param', candidate_params(c,:), 'name', candidate_names{c});

    fprintf('  %s: TF_peak=%.4f, stable=%.1f%%\n', candidate_names{c}, ...
        max(TF_c_lin, [], 'omitnan'), ...
        100 * sum(maxMu_c < tol_stable) / max(1, N_pts));
end

%% --- 4.1 Wang BG 模型基线（纯机械，lam=0, same mechanical params）---
fprintf('\n--- Wang BG 模型基线 (K1=1, K2=0, lam=0) ---\n');
sysP_wang = sysP0;
sysP_wang(7)  = 0.0;   % lam = 0 (no circuit)
sysP_wang(8)  = 0.0;   % kap_e = 0
sysP_wang(9)  = 0.0;   % kap_c = 0
sysP_wang(10) = 0.0;   % sigma = 0

try
    [Om_wang, TF_dB_wang, x_res_wang] = call_robust_frf(sysP_wang, Fw_opt, 5000);
catch
    Om_wang = []; TF_dB_wang = []; x_res_wang = [];
end

if ~isempty(Om_wang) && length(Om_wang) >= 10
    TF_wang_lin = 10.^(TF_dB_wang/20);
    ok_wang = true(size(Om_wang));

    % Floquet 扫描 Wang 基线
    N_floq_wang = min(200, length(Om_wang));
    floq_idx_w = unique(round(linspace(1, length(Om_wang), N_floq_wang)));
    maxMu_wang = nan(N_floq_wang, 1);
    Om_floq_wang = Om_wang(floq_idx_w);

    for k = 1:N_floq_wang
        xc_k = x_res_wang(1:15, floq_idx_w(k));
        Om_k = x_res_wang(16, floq_idx_w(k));
        maxMu_wang(k) = compute_floquet_fast(xc_k, sysP_wang, Om_k, Nt_floquet);
    end

    % 添加 Wang 基线为第 4 个候选
    results{4} = struct('TF', TF_wang_lin, 'TF_dB', TF_dB_wang, ...
        'maxMu', maxMu_wang, 'ok', ok_wang, ...
        'Om', Om_wang, 'Om_floq', Om_floq_wang, ...
        'x_res', x_res_wang, ...
        'param', [0, 0, 0], 'name', 'Wang BG baseline (K1=1,K2=0,lam=0)');
    x_res_all{4} = x_res_wang;
    Nc_total = 4;
    fprintf('  Wang BG baseline: TF_peak=%.4f (%.1f dB), stable=%.1f%%\n', ...
        max(TF_wang_lin,[],'omitnan'), 20*log10(max(TF_wang_lin,[],'omitnan')), ...
        100 * sum(maxMu_wang < tol_stable) / max(1, N_floq_wang));
else
    Nc_total = 3;
    fprintf('  Wang BG baseline: arc-length failed, 仅 3 候选对比\n');
end

%% --- 4.2 选择综合最优候选 ---
best_for_plot = 3;  % 默认回退到 Ph2 raw
best_tf_for_sel = inf;
for c_sel = [1, 2, 3]  % 优先检查 Ph3 候选
    maxMu_c = results{c_sel}.maxMu;
    if isempty(maxMu_c) || all(isnan(maxMu_c))
        continue;
    end
    stable_pct = 100 * sum(maxMu_c < tol_stable) / max(1, length(maxMu_c));
    TF_c = max(results{c_sel}.TF, [], 'omitnan');
    if stable_pct > 90 && TF_c < best_tf_for_sel
        best_tf_for_sel = TF_c;
        best_for_plot = c_sel;
    end
end

fprintf('\n  选用 %s 绘图 (综合最优)\n', candidate_names{best_for_plot});

sysP_best = sysP0;
sysP_best(8)  = candidate_params(best_for_plot, 2);
sysP_best(9)  = candidate_params(best_for_plot, 3);
sysP_best(10) = candidate_params(best_for_plot, 1);
TF_best   = results{best_for_plot}.TF;
TF_dB_best = results{best_for_plot}.TF_dB;
Om_best   = results{best_for_plot}.Om;
maxMu_best = results{best_for_plot}.maxMu;
Om_floq_best = results{best_for_plot}.Om_floq;
ok_best    = results{best_for_plot}.ok;
best_p_final = candidate_params(best_for_plot, :);

%% --- 4.3 纯机械基线 (lam=0, no circuit) ---
sysP_base = sysP0;
sysP_base(7)  = 0.0;   % lam = 0
sysP_base(8)  = 0.0;   % kap_e = 0
sysP_base(9)  = 0.0;   % kap_c = 0
sysP_base(10) = 0.0;   % sigma = 0

try
    [Om_base, TF_dB_base] = arc_length_frf(sysP_base, 10.0, ...
        'Fw', Fw_opt, 'Step', -0.01, 'Steps', 3000);
catch
    Om_base = []; TF_dB_base = [];
end

if isempty(Om_base)
    TF_base = nan; ok_base = false; Om_base = nan;
else
    TF_base = 10.^(TF_dB_base/20);
    ok_base = true(size(Om_base));
end

fprintf('\n');
[TFpk_b, ipk] = max(TF_best);
[TFpk_base, ipk_base] = max(TF_base);
TFpk_b_dB = 20*log10(TFpk_b);
TFpk_base_dB = 20*log10(TFpk_base);

%% --- 4.4 Wang 2017 四项性能指标计算 ---
fprintf('\n--- Wang 2017 性能指标计算 ---\n');
wang_data = cell(Nc_total, 1);
wang_names = cell(Nc_total, 1);

for c = 1:Nc_total
    wang_names{c} = results{c}.name;
    if isfield(results{c}, 'x_res') && ~isempty(results{c}.x_res)
        xr = results{c}.x_res;
        Om_c = results{c}.Om;
        TF_dB_c = results{c}.TF_dB;

        % 提取有效点
        valid_c = isfinite(Om_c) & isfinite(TF_dB_c) & (Om_c > 0);
        Om_v = Om_c(valid_c);
        TF_v = TF_dB_c(valid_c);

        % 排序
        [Om_s, srt] = sort(Om_v);
        TF_s = TF_v(srt);

        % I1: 峰值动态位移
        x1_h = xr(1:5, :)';
        amp1 = sqrt(x1_h(:,2).^2 + x1_h(:,3).^2);
        amp1 = amp1(valid_c);
        I1 = max(amp1(srt));

        % I2: 峰值力传递率 (dB)
        I2 = max(TF_s);

        % I3: 0 dB 穿越频率
        I3 = find_cross_frequency(Om_s, TF_s, 0);

        % I4: -40 dB 穿越频率
        I4 = find_cross_frequency(Om_s, TF_s, -40);

        wang_data{c} = struct('I1', I1, 'I2', I2, 'I3', I3, 'I4', I4);
        fprintf('  %-35s: I1=%.4f, I2=%.2f dB, I3=%s, I4=%s\n', ...
            wang_names{c}, I1, I2, ...
            iff(isnan(I3), 'N/A', sprintf('%.4f', I3)), ...
            iff(isnan(I4), 'N/A', sprintf('%.4f', I4)));
    else
        wang_data{c} = struct('I1', nan, 'I2', nan, 'I3', nan, 'I4', nan);
    end
end

%% --- 4.5 算子 K(Ω) 分解 (解释最优参数为何有效) ---
fprintf('\n--- 算子解释：最优参数的物理机制 ---\n');
Om_op = logspace(log10(0.05), log10(10), 500).';
den_k = best_p_final(2)*Om_op.^2 - 1i*best_p_final(1)*Om_op - best_p_final(3);
K_opt = (theta^2) .* Om_op.^2 ./ den_k;
Kr_opt = real(K_opt);
Ki_opt = imag(K_opt);
Ceq_opt = Ki_opt ./ max(Om_op, 1e-12);
Meq_opt = (theta^2)/best_p_final(3) * ones(size(Om_op));
Keq_opt = Kr_opt + Om_op.^2 .* Meq_opt;
Om_e = sqrt(abs(best_p_final(3)/best_p_final(2)));  % 电路电气谐振频率

% 算子物理解释
sigma_opt = best_p_final(1);
kap_c_opt = best_p_final(3);
mech_damping = 2 * ze1 * sqrt(be1 + al1);  % 机械阻尼近似
if sigma_opt < 0
    fprintf('  sigma=%.4f < 0 → NIC 注入"负阻尼", 部分抵消机械阻尼\n', sigma_opt);
    fprintf('    机械阻尼 ~%.4f, NIC 有效负阻 = %.4f\n', mech_damping, abs(sigma_opt));
end
if kap_c_opt < 0
    fprintf('  kap_c=%.4f < 0 → Meq = theta^2/kap_c = %.4f < 0 (负等效惯容)\n', ...
        kap_c_opt, Meq_opt(1));
    fprintf('    负 Meq 降低 Keq 在低频的值, 助于拓宽隔离带\n');
end
fprintf('  电路特征频率 Omega_e = %.4f (电气谐振)\n', Om_e);

%% --- 4.6 图 1: 力传递率 + Floquet 对比 (4-panel) ---
fontName = 'Times New Roman';
fsLab = 12; fsTit = 13;

figure('Color','w','Position',[30 30 1400 900]);

% Panel A: TF 对比 (dB)
subplot(2,3,1);
semilogx(Om_best, TF_dB_best, 'b-', 'LineWidth', 1.8); hold on;
if ~all(isnan(TF_base))
    semilogx(Om_base, TF_dB_base, 'Color', [0.6 0.6 0.6], ...
        'LineWidth', 1.5, 'LineStyle', '--');
end
if exist('TF_dB_wang', 'var') && ~isempty(TF_dB_wang)
    semilogx(Om_wang, TF_dB_wang, 'r-', 'LineWidth', 1.5);
end
grid on; box on;
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('T_F (dB)', 'FontName', fontName, 'FontSize', fsLab);
title('(A) Force Transmissibility (dB)', 'FontName', fontName, 'FontSize', fsTit);
yline(0, 'k--', 'LineWidth', 1.0);
leg_strs = {'Optimized EMSD', 'Pure mechanical (lam=0)'};
if exist('TF_dB_wang', 'var') && ~isempty(TF_dB_wang)
    leg_strs{3} = 'Wang BG (K1=1,K2=0)';
end
legend(leg_strs, 'Location', 'best', 'FontName', fontName, 'FontSize', 9);

reduction = (TFpk_base - TFpk_b) / TFpk_base * 100;
text(Om_best(ipk), TF_dB_best(ipk), ...
    sprintf('  Peak: %.1f dB\n  \\Delta: %.1f%%', TF_dB_best(ipk), reduction), ...
    'FontName', fontName, 'FontSize', 9);

% Panel B: TF 对比 (线性)
subplot(2,3,2);
semilogx(Om_best, TF_best, 'b-', 'LineWidth', 1.8); hold on;
if ~all(isnan(TF_base))
    semilogx(Om_base, TF_base, 'Color', [0.6 0.6 0.6], ...
        'LineWidth', 1.5, 'LineStyle', '--');
end
if exist('TF_wang_lin', 'var') && ~isempty(TF_wang_lin)
    semilogx(Om_wang, TF_wang_lin, 'r-', 'LineWidth', 1.5);
end
grid on; box on;
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('T_F (linear)', 'FontName', fontName, 'FontSize', fsLab);
title('(B) Force Transmissibility (linear)', 'FontName', fontName, 'FontSize', fsTit);
legend(leg_strs, 'Location', 'best', 'FontName', fontName, 'FontSize', 9);

% Panel C: Floquet 最优候选
subplot(2,3,3);
semilogx(Om_floq_best, maxMu_best, 'r.-', 'LineWidth', 1.2, 'MarkerSize', 6); hold on;
yline(tol_stable, 'k--', 'LineWidth', 1.2);
yline(1.0, ':', 'Color', [0.4 0.4 0.4]);
grid on; box on;
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('max|\mu|', 'FontName', fontName, 'FontSize', fsLab);
title(sprintf('(C) Floquet: %s', candidate_names{best_for_plot}), ...
    'FontName', fontName, 'FontSize', fsTit);
ylim([0, max(1.6, 1.1*max(maxMu_best))]);

% Panel D: Displacement amplitude 对比
subplot(2,3,4);
xr_best = results{best_for_plot}.x_res;
x1h_best = xr_best(1:5, :)';
amp_best = sqrt(x1h_best(:,2).^2 + x1h_best(:,3).^2);
valid_b = isfinite(Om_best) & isfinite(amp_best) & (Om_best > 0);
loglog(Om_best(valid_b), amp_best(valid_b), 'b-', 'LineWidth', 1.8); hold on;

if exist('x_res_wang', 'var') && ~isempty(x_res_wang)
    x1h_w = x_res_wang(1:5, :)';
    amp_w = sqrt(x1h_w(:,2).^2 + x1h_w(:,3).^2);
    valid_w = isfinite(Om_wang) & isfinite(amp_w) & (Om_wang > 0);
    loglog(Om_wang(valid_w), amp_w(valid_w), 'r-', 'LineWidth', 1.5);
end
grid on; box on;
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('|X_1| (fund.)', 'FontName', fontName, 'FontSize', fsLab);
title('(D) Displacement Amplitude', 'FontName', fontName, 'FontSize', fsTit);

% Panel E: Floquet 全候选对比
subplot(2,3,5);
for c = 1:Nc_total
    if isfield(results{c}, 'Om_floq') && ~all(isnan(results{c}.maxMu))
        semilogx(results{c}.Om_floq, results{c}.maxMu, ...
            'LineWidth', 1.2, 'DisplayName', wang_names{c}); hold on;
    end
end
yline(tol_stable, 'k--', 'LineWidth', 1.3);
yline(1.0, ':', 'Color', [0.4 0.4 0.4]);
grid on; box on;
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('max|\mu|', 'FontName', fontName, 'FontSize', fsLab);
title('(E) Floquet Comparison', 'FontName', fontName, 'FontSize', fsTit);
ylim([0, 2.0]);
legend('Location', 'best', 'FontName', fontName, 'FontSize', 8);

% Panel F: K operator complex plane
subplot(2,3,6);
scatter(Kr_opt, Ki_opt, 15, Om_op, 'filled'); hold on;
colormap(jet); cb = colorbar; cb.Label.String = '\Omega';
grid on; box on;
xlabel('K_r (real)', 'FontName', fontName, 'FontSize', fsLab);
ylabel('K_i (imag)', 'FontName', fontName, 'FontSize', fsLab);
title(sprintf('(F) K(\\Omega) Complex Plane\n\\sigma=%.3f, \\kappa_e=%.3f, \\kappa_c=%.3f', ...
    best_p_final(1), best_p_final(2), best_p_final(3)), ...
    'FontName', fontName, 'FontSize', fsTit);
axis equal;

sgtitle(sprintf(['Optimization Results: Best = \\sigma=%.4f, \\kappa_e=%.4f, ' ...
    '\\kappa_c=%.4f | Peak TF=%.4f (%.1f%% reduction)'], ...
    best_p_final(1), best_p_final(2), best_p_final(3), TFpk_b, reduction), ...
    'FontName', fontName, 'FontSize', 14);

%% --- 4.7 图 2: 算子解释 (6-panel) ---
figure('Color','w','Position',[50 50 1400 900]);

% (A) K complex plane (Kr vs Ki, colored by Omega)
subplot(3,2,1);
scatter(Kr_opt, Ki_opt, 12, Om_op, 'filled'); hold on;
colormap(jet); cb = colorbar; cb.Label.String = '\Omega';
xlabel('K_r', 'FontName', fontName, 'FontSize', fsLab);
ylabel('K_i', 'FontName', fontName, 'FontSize', fsLab);
title('(A) K(\Omega) Complex Plane', 'FontName', fontName, 'FontSize', fsTit);
grid on; box on;

% (B) Ceq = imag(K)/Omega
subplot(3,2,2);
semilogx(Om_op, Ceq_opt, 'LineWidth', 1.8);
grid on; box on;
xline(Om_e, 'k--', 'LineWidth', 1.0);
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('C_{eq}(\Omega)', 'FontName', fontName, 'FontSize', fsLab);
title('(B) Equivalent Damping C_{eq}', 'FontName', fontName, 'FontSize', fsTit);

% (C) Keq(Omega)
subplot(3,2,3);
semilogx(Om_op, Keq_opt, 'LineWidth', 1.5);
grid on; box on;
xline(Om_e, 'k--', 'LineWidth', 1.0);
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('K_{eq}(\Omega)', 'FontName', fontName, 'FontSize', fsLab);
title('(C) Equivalent Stiffness K_{eq}', 'FontName', fontName, 'FontSize', fsTit);

% (D) Meq = theta^2/kap_c
subplot(3,2,4);
semilogx(Om_op, Meq_opt, 'LineWidth', 1.5);
grid on; box on;
xline(Om_e, 'k--', 'LineWidth', 1.0);
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('M_{eq}(\Omega)', 'FontName', fontName, 'FontSize', fsLab);
title(sprintf('(D) Equivalent Inertia M_{eq}=%.4f', Meq_opt(1)), ...
    'FontName', fontName, 'FontSize', fsTit);

% (E) Force transmissibility (dB) for best candidate
subplot(3,2,5);
semilogx(Om_best, TF_dB_best, 'b-', 'LineWidth', 1.8); hold on;
if ~all(isnan(TF_base))
    semilogx(Om_base, TF_dB_base, 'Color', [0.6 0.6 0.6], ...
        'LineWidth', 1.5, 'LineStyle', '--');
end
grid on; box on;
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('T_F (dB)', 'FontName', fontName, 'FontSize', fsLab);
title('(E) Force Transmissibility', 'FontName', fontName, 'FontSize', fsTit);
yline(0, 'k--', 'LineWidth', 1.0);
legend({'Optimal', 'No circuit'}, 'Location', 'best', 'FontName', fontName, 'FontSize', 9);

% (F) Displacement amplitude
subplot(3,2,6);
loglog(Om_best(valid_b), amp_best(valid_b), 'b-', 'LineWidth', 1.8);
grid on; box on;
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('|X_1|', 'FontName', fontName, 'FontSize', fsLab);
title('(F) Displacement Amplitude', 'FontName', fontName, 'FontSize', fsTit);

sgtitle(sprintf(['Operator Interpretation: \\sigma=%.4f, \\kappa_e=%.4f, \\kappa_c=%.4f | ' ...
    '\\Omega_e=%.3f, M_{eq}=%.3f'], ...
    best_p_final(1), best_p_final(2), best_p_final(3), Om_e, Meq_opt(1)), ...
    'FontName', fontName, 'FontSize', 14);

%% --- 4.8 图 3: Wang 风格参数扫掠 (3×2) ---
fprintf('\n--- Wang 风格参数扫掠 ---\n');
N_swp = 15;  % 每个参数的扫掠点数
Budget_swp = 2000;

% 扫掠范围: 以最优值为中心 ± 边界
sigma_swp = linspace(max(lb(1), best_p_final(1)-2.5), min(ub(1), best_p_final(1)+2.5), N_swp)';
kap_e_swp = linspace(max(lb(2), best_p_final(2)*0.1),  min(ub(2), best_p_final(2)*3), N_swp)';
kap_c_swp = linspace(max(lb(3), best_p_final(3)-2.5), min(ub(3), best_p_final(3)+2.5), N_swp)';

% 初始化扫掠结果
swp_I1 = cell(3, N_swp); swp_I2 = cell(3, N_swp);
swp_I3 = cell(3, N_swp); swp_I4 = cell(3, N_swp);

% 扫掠 sigma
fprintf('  扫掠 sigma ...\n');
for j = 1:N_swp
    p_s = best_p_final; p_s(1) = sigma_swp(j);
    if abs(p_s(3)) < 1e-4, continue; end
    try
        [swp_I1{1,j}, swp_I2{1,j}, swp_I3{1,j}, swp_I4{1,j}] = ...
            eval_one_param(p_s, sysP0, Fw_opt, Budget_swp, Nt_floquet);
    catch
    end
end

% 扫掠 kap_e
fprintf('  扫掠 kap_e ...\n');
for j = 1:N_swp
    p_s = best_p_final; p_s(2) = kap_e_swp(j);
    if abs(p_s(3)) < 1e-4, continue; end
    try
        [swp_I1{2,j}, swp_I2{2,j}, swp_I3{2,j}, swp_I4{2,j}] = ...
            eval_one_param(p_s, sysP0, Fw_opt, Budget_swp, Nt_floquet);
    catch
    end
end

% 扫掠 kap_c
fprintf('  扫掠 kap_c ...\n');
for j = 1:N_swp
    p_s = best_p_final; p_s(3) = kap_c_swp(j);
    if abs(p_s(3)) < 1e-4, continue; end
    try
        [swp_I1{3,j}, swp_I2{3,j}, swp_I3{3,j}, swp_I4{3,j}] = ...
            eval_one_param(p_s, sysP0, Fw_opt, Budget_swp, Nt_floquet);
    catch
    end
end

fprintf('  参数扫掠完成\n');

% 绘制图 3
figure('Color','w','Position',[80 80 1100 900]);
swp_params = {sigma_swp, kap_e_swp, kap_c_swp};
swp_labels = {'\sigma (resistance ratio)', '\kappa_e (inductance)', '\kappa_c (capacitance^-1)'};
subplot_titles = {'I1: Peak Displacement', 'I2: Peak TF (dB)', ...
                  'I3: 0dB Cross Frequency', 'I4: -40dB Cross Frequency'};

% Row 1: sigma sweep → col 1: I1+I2, col 2: I3+I4
% Row 2: kap_e sweep → same
% Row 3: kap_c sweep → same
for row = 1:3
    par_vals = swp_params{row};

    % Col 1: I1 + I2
    subplot(3, 2, 2*(row-1)+1);
    yyaxis left;
    I1_row = safe_get_wang_field(swp_I1(row,:), 'I1');
    plot(par_vals, I1_row, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6);
    ylabel('I1: Peak Disp.', 'FontName', fontName, 'FontSize', fsLab-2);
    yyaxis right;
    I2_row = safe_get_wang_field(swp_I2(row,:), 'I2');
    plot(par_vals, I2_row, 'rs-', 'LineWidth', 1.5, 'MarkerSize', 6);
    ylabel('I2: Peak TF (dB)', 'FontName', fontName, 'FontSize', fsLab-2);
    grid on; box on;
    xlabel(swp_labels{row}, 'FontName', fontName, 'FontSize', fsLab);
    title(sprintf('%s → I1 + I2', swp_labels{row}), ...
        'FontName', fontName, 'FontSize', fsTit-1);
    xline(best_p_final(row), 'k--', 'LineWidth', 1.0);

    % Col 2: I3 + I4
    subplot(3, 2, 2*(row-1)+2);
    yyaxis left;
    I3_row = safe_get_wang_field(swp_I3(row,:), 'I3');
    plot(par_vals, I3_row, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6);
    ylabel('I3: 0dB Cross \Omega', 'FontName', fontName, 'FontSize', fsLab-2);
    yyaxis right;
    I4_row = safe_get_wang_field(swp_I4(row,:), 'I4');
    plot(par_vals, I4_row, 'rs-', 'LineWidth', 1.5, 'MarkerSize', 6);
    ylabel('I4: -40dB Cross \Omega', 'FontName', fontName, 'FontSize', fsLab-2);
    grid on; box on;
    xlabel(swp_labels{row}, 'FontName', fontName, 'FontSize', fsLab);
    title(sprintf('%s → I3 + I4', swp_labels{row}), ...
        'FontName', fontName, 'FontSize', fsTit-1);
    xline(best_p_final(row), 'k--', 'LineWidth', 1.0);
end

sgtitle(sprintf(['Parameter Sweep Around Optimum: \\sigma=%.3f, \\kappa_e=%.3f, \\kappa_c=%.3f'], ...
    best_p_final(1), best_p_final(2), best_p_final(3)), ...
    'FontName', fontName, 'FontSize', 14);

%% --- 4.9 Wang 2017 性能指标对比表 ---
fprintf('\n============= Wang 2017 性能指标对比表 =============\n');
fprintf('%-40s %-16s %-14s %-16s %-16s\n', ...
    'Candidate', 'I1:PeakDisp', 'I2:PeakTF', 'I3:0dB Cross', 'I4:-40dB Cross');
fprintf('%s\n', repmat('-', 1, 104));
for c = 1:Nc_total
    wd = wang_data{c};
    I3_str = iff(isnan(wd.I3), 'N/A', sprintf('%.4f', wd.I3));
    I4_str = iff(isnan(wd.I4), 'N/A', sprintf('%.4f', wd.I4));
    fprintf('%-40s %-16.4f %-14.2f %-16s %-16s\n', ...
        wang_names{c}, wd.I1, wd.I2, I3_str, I4_str);
end
fprintf('==========================================================\n');

%% ========================================================================
% 5. 结果汇总
%% ========================================================================
% 数据一致性单一来源 (Single Source of Truth):
%   - 基线定义: 纯机械系统 (lambda=0, 所有电路参数为零) 在相同 Fw=0.005,
%               Omega 范围 [0.05, 10.0] 内的 TF 峰值（弧长延拓）。
%   - 峰值降低百分比 = (baseline_TF_peak - opt_TF_peak) / baseline_TF_peak * 100%
%   - 此数值必须在以下位置保持一致：摘要、引言、第5章结果、第7章结论。
%   - 在修改任何模型参数后，请重新运行本脚本并更新下方注释。
fprintf('\n==============================================\n');
fprintf('           最终优化结果汇总\n');
fprintf('==============================================\n');
fprintf('  选用方案: %s\n', candidate_names{best_for_plot});
fprintf('  sigma (电阻比)     = %.6f\n', best_p_final(1));
fprintf('  kap_e (电感)       = %.6f\n', best_p_final(2));
fprintf('  kap_c (电容倒数)   = %.6f\n', best_p_final(3));
fprintf('  Omega_e (电路特征) = %.6f\n', Om_e);
fprintf('  TF_peak (最优)     = %.6f (%.2f dB)\n', TFpk_b, TFpk_b_dB);
fprintf('  TF_peak (baseline) = %.6f (%.2f dB)\n', TFpk_base, TFpk_base_dB);
fprintf('  峰值降低            = %.1f%%\n', reduction);
fprintf('  Meq (等效惯容)     = %.6f\n', Meq_opt(1));
fprintf('  Ceq 峰值           = %.6f\n', max(abs(Ceq_opt)));
fprintf('  稳定点比例          = %.1f%%\n', ...
    100 * sum(maxMu_best < tol_stable) / max(1, length(maxMu_best)));

% --- 5.1 NIC 有源功率评估 (Supplement 2) ---
% Compute NIC active power at design excitation
sigma_active = 1.0 - best_p_final(1);
% Find the TF peak frequency for power computation
[~, idx_peak] = max(TF_best);
Om_peak = Om_best(idx_peak);
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

% --- 5.2 全候选对比 ---
fprintf('  全候选对比:\n');
for c = 1:Nc_total
    TF_c = results{c}.TF; maxMu_c = results{c}.maxMu;
    param_c = results{c}.param;
    fprintf('    %s: sigma=%.4f kap_e=%.4f kap_c=%.4f  TFpk=%.4f  stb=%.1f%%\n', ...
        wang_names{c}, param_c(1), param_c(2), ...
        param_c(3), max(TF_c,[],'omitnan'), ...
        100*sum(maxMu_c<tol_stable)/max(1,length(maxMu_c)));
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
    'candidate_names', 'candidate_params', 'results', 'wang_data', 'wang_names', ...
    'TFpk_b', 'TFpk_base', 'reduction', 'Om_e', 'Meq_opt', ...
    'Ceq_opt', 'Keq_opt', 'Kr_opt', 'Ki_opt', ...
    'sysP0', 'Fw_opt', 'Om_min', 'Om_max', 'tol_stable', ...
    'lb', 'ub', 'Nc_total');

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

function [Om, TF_dB, x_res] = call_robust_frf(sysP, Fw_val, budget)
%% Robust FRF evaluation: upward sweep primary, downward sweep fallback
    try
        [Om, TF_dB, x_res] = arc_length_frf_robust(sysP, 10.0, ...
            'Fw', Fw_val, 'Budget', budget, 'Verbose', false);
        if ~isempty(Om) && length(Om) >= 20
            return;
        end
    catch
    end
    % Fallback: standard downward sweep (suppress verbose output with evalc)
    try
        steps = max(budget, 1500);
        [~, Om, TF_dB, x_res] = evalc( ...
            'arc_length_frf(sysP, 10.0, ''Fw'', Fw_val, ''Step'', -0.01, ''Steps'', steps)');
    catch
        Om = []; TF_dB = []; x_res = [];
    end
end

function [J, out] = objective_wrapper(z, lb, ub, sysP0, Fw_val, ...
                                      Nt_flo, tol_stab)
% Phase 3 (fminsearch) 目标函数 — 使用弧长延拓 FRF 评估
    p = lb + (ub - lb) .* sigmoid(z);
    sigma_i = p(1); kap_e_i = p(2); kap_c_i = p(3);

    sysP = sysP0;
    sysP(8)  = kap_e_i;
    sysP(9)  = kap_c_i;
    sysP(10) = sigma_i;

    % --- 弧长延拓 FRF ---
    global ParamMin ParamMax
    ParamMin = 0.05;
    ParamMax = 10.5;

    try
        [Om, TF_dB, x_res] = call_robust_frf(sysP, Fw_val, 2000);
    catch
        J = 1e8;
        out = struct('TF_peak', inf);
        return;
    end

    if isempty(Om) || isempty(TF_dB)
        J = 1e6;
        out = struct('TF_peak', inf);
        return;
    end

    % 线性 TF 峰值
    TF_lin = 10.^(TF_dB/20);
    [peakTF, idx_peak] = max(TF_lin);

    % --- Floquet 稳定性抽查（低频/峰值/高频 3 个去重点）---
    Om_raw = x_res(16, :)';
    valid_x = Om_raw > 0 & isfinite(Om_raw);
    x_res_f = x_res(:, valid_x);
    npts = size(x_res_f, 2);

    if npts >= 3
        Om_peak_val = Om(idx_peak);
        [~, idx_peak_x] = min(abs(x_res_f(16, :)' - Om_peak_val));
        check_idx = unique([1, idx_peak_x, npts]);
        mu_vals = zeros(1, length(check_idx));
        for k = 1:length(check_idx)
            xc_k = x_res_f(1:15, check_idx(k));
            Om_k = x_res_f(16, check_idx(k));
            mu_vals(k) = compute_floquet_fast(xc_k, sysP, Om_k, Nt_flo);
        end
        maxMu = max(mu_vals);

        % 稳定性惩罚
        pen_stab = 0;
        if maxMu > tol_stab
            pen_stab = 5e1 * (maxMu - tol_stab)^2 + 1e1;
        end
    else
        pen_stab = 0;
    end

    J = peakTF + pen_stab;
    out = struct('TF_peak', peakTF);
end

function s = sigmoid(z)
    s = 1 ./ (1 + exp(-z));
end

function z = inv_sigmoid(s)
    s = min(max(s, 1e-6), 1 - 1e-6);
    z = log(s ./ (1 - s));
end

function cb = make_early_stop(max_stall)
% 创建 fminsearch 提前终止回调
% 如果连续 max_stall 次迭代目标函数未改善 >1%，则停止
    best_fval = inf;
    stall_count = 0;

    cb = @early_stop_callback;

    function stop = early_stop_callback(~, optimValues, state)
        switch state
            case 'init'
                best_fval = optimValues.fval;
                stall_count = 0;
            case 'iter'
                if optimValues.fval < best_fval * 0.99  % 1% 改善
                    best_fval = optimValues.fval;
                    stall_count = 0;
                else
                    stall_count = stall_count + 1;
                end
        end
        stop = (stall_count >= max_stall);
    end
end

%% --- 新增辅助函数 ---

function [I1, I2, I3, I4] = eval_one_param(params, sysP0, Fw_val, budget, Nt_floquet)
%% 评估单组参数的 Wang 四项指标
% Returns: I1 (peak displacement), I2 (peak TF dB), I3 (0dB cross), I4 (-40dB cross)
% Returns empty if evaluation fails
    sysP = sysP0;
    sysP(8)  = params(2);
    sysP(9)  = params(3);
    sysP(10) = params(1);

    try
        [Om, TF_dB, x_res] = call_robust_frf(sysP, Fw_val, budget);
    catch
        I1 = []; I2 = []; I3 = []; I4 = []; return;
    end

    if isempty(Om) || length(Om) < 10
        I1 = []; I2 = []; I3 = []; I4 = []; return;
    end

    % 筛选有效点
    valid = isfinite(Om) & isfinite(TF_dB) & (Om > 0);
    Om_v = Om(valid);
    TF_v = TF_dB(valid);

    % 排序
    [Om_s, srt] = sort(Om_v);
    TF_s = TF_v(srt);

    % I1: 峰值位移
    x1_h = x_res(1:5, :)';
    amp1 = sqrt(x1_h(:,2).^2 + x1_h(:,3).^2);
    amp1 = amp1(valid);
    I1 = max(amp1(srt));

    % I2: 峰值 TF (dB)
    I2 = max(TF_s);

    % I3: 0 dB 穿越
    I3 = find_cross_frequency(Om_s, TF_s, 0);

    % I4: -40 dB 穿越
    I4 = find_cross_frequency(Om_s, TF_s, -40);
end

function f_cross = find_cross_frequency(f, y, level)
%% 寻找第一次从上方穿越给定水平线的频率
% Input:
%   f     : 频率升序数组
%   y     : 对应的传递率 dB 值
%   level : 目标 dB 值
% Output:
%   f_cross : 穿越频率，若无穿越则返回 NaN
    % 找到第一个 <= level 的点
    idx = find(y <= level, 1, 'first');
    if isempty(idx)
        f_cross = NaN;
        return;
    end

    if idx == 1
        f_cross = f(1);
    else
        f1 = f(idx-1); y1 = y(idx-1);
        f2 = f(idx);   y2 = y(idx);
        if y1 == y2
            f_cross = f1;
        else
            f_cross = f1 + (level - y1) * (f2 - f1) / (y2 - y1);
        end
    end
end

function [Om_low, Om_high] = find_isolation_band(Om, TF_dB, threshold)
%% 寻找 TF_dB 连续低于 threshold 的隔离频带
% Input:
%   Om        : 升序频率数组
%   TF_dB     : 对应 TF dB 值
%   threshold : 隔离阈值 (dB)
% Output:
%   Om_low  : 隔离带起始频率 (first crossing below threshold)
%   Om_high : 隔离带结束频率 (last crossing below threshold before exceeding)
%   If no isolation band, both are NaN.
    idx_below = find(TF_dB <= threshold);
    if isempty(idx_below)
        Om_low = NaN; Om_high = NaN;
        return;
    end
    Om_low = Om(idx_below(1));
    % Find longest continuous below-threshold stretch
    diffs = diff(idx_below);
    breaks = find(diffs > 1);
    if isempty(breaks)
        Om_high = Om(idx_below(end));
    else
        % Take the longest segment
        segments = [idx_below(1); idx_below(breaks+1)];
        ends = [idx_below(breaks); idx_below(end)];
        lens = ends - segments + 1;
        [~, best] = max(lens);
        Om_low  = Om(segments(best));
        Om_high = Om(ends(best));
    end
end

function result = safe_get_wang_field(cell_row, field_name)
%% Robustly extract a field from cell array elements, returning nan on failure
    result = nan(1, length(cell_row));
    for j = 1:length(cell_row)
        x = cell_row{j};
        if isempty(x) || ~isstruct(x) || ~isfield(x, field_name)
            result(j) = nan;
        else
            result(j) = x.(field_name);
        end
    end
end

function v = iff(cond, a, b)
%% Inline if: v = cond ? a : b
    if cond, v = a; else, v = b; end
end
