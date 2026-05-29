%% =========================================================
% verify_peak2_by_changing_Csh.m
%
% 目的：
%   验证第二个峰值附近，电路固有频率 Omega_e 的匹配效果。
%
% 核心思想：
%   机电耦合系数 lambda 不人为调节，而由物理参数确定：
%
%       lambda = Kt*Ke*wn/(k1*R0)
%
%   固定 VCM、机械结构、线圈内阻、电阻、电感：
%
%       Kt, Ke, m1, k1, R0, Rt, Lsh 固定
%
%   只改变外接电容 Csh，使电路固有频率 Omega_e
%   移动到 baseline 第二峰附近。
%
%   Omega_e = 1/(wn*sqrt(Lsh*Csh))
%
%   因此：
%
%       Csh = 1/(Lsh*wn^2*Omega_e^2)
%
% 依赖：
%   nondim_temp2.m
%   newton.m
%   branch_follow2.m
%   branch_aux2.m
%
% 注意：
%   你当前 nondim_temp2.m 中，当 lam≈0 时会自动加入层间阻尼 zeta12=0.05；
%   当 lam>0 时该层间阻尼关闭。
%   因此 baseline 与接电路曲线的机械阻尼条件不完全一致。
%   正式论文中建议把 zeta12 独立写入 sysP。
% =========================================================

clc; clear; close all;

%% =========================================================
% 0. 全局变量
%% =========================================================
global Fw FixedOmega ParamMin ParamMax

Fw = 0.005;
FixedOmega = [];     % 扫频模式：nondim_temp2 中 y(16)=Omega

%% =========================================================
% 1. 物理电路与 VCM 参数
%% =========================================================
Kt = 7.474;
Ke = 7.474;
m1 = 2.2;
k1 = 3000;
R0 = 3.8;

wn = sqrt(k1/m1);

Rt  = 2.3674;
Lsh = 0.04065;

% 由物理参数确定的无量纲电路参数
lambda_phys = Kt*Ke*wn/(k1*R0);
kap_e_phys  = Lsh*wn/R0;
sigma_phys  = Rt/R0;

fprintf('\n===== Fixed physical electromechanical parameters =====\n');
fprintf('Kt       = %.6f\n', Kt);
fprintf('Ke       = %.6f\n', Ke);
fprintf('m1       = %.6f kg\n', m1);
fprintf('k1       = %.6f N/m\n', k1);
fprintf('R0       = %.6f Ohm\n', R0);
fprintf('wn       = %.6f rad/s\n', wn);
fprintf('Rt       = %.6f Ohm\n', Rt);
fprintf('Lsh      = %.6f H\n', Lsh);
fprintf('lambda   = %.8f\n', lambda_phys);
fprintf('kap_e    = %.8f\n', kap_e_phys);
fprintf('sigma    = %.8f\n', sigma_phys);

%% =========================================================
% 2. 机械 / QZS 参数：与你当前设置保持一致
%% =========================================================
mu   = 0.2;     % 质量比 m2/m1
beta = 2.0;     % 下层竖向线性刚度比
K1   = 1.0;     % 上层水平弹簧刚度比
K2   = 0.5;     % 下层水平弹簧刚度比
U    = 2.0;     % 几何非线性尺度参数
L    = 4/9;     % QZS 长度比

% 由 L=4/9, K1=1, alpha1=0 反推
v = 2.5;

% Wang 写法：
% alpha = v - 2*K*(1-L)/L
alpha1 = v    - 2*K1*(1-L)/L;
alpha2 = beta - 2*K2*(1-L)/L;

gamma1 = K1/(U^2 * L^3);
gamma2 = K2/(U^2 * L^3);

fprintf('\n===== Mechanical / QZS parameters =====\n');
fprintf('mu       = %.6f\n', mu);
fprintf('beta     = %.6f\n', beta);
fprintf('K1       = %.6f\n', K1);
fprintf('K2       = %.6f\n', K2);
fprintf('U        = %.6f\n', U);
fprintf('L        = %.6f\n', L);
fprintf('v        = %.6f\n', v);
fprintf('alpha1   = %.8f\n', alpha1);
fprintf('alpha2   = %.8f\n', alpha2);
fprintf('gamma1   = %.8f\n', gamma1);
fprintf('gamma2   = %.8f\n', gamma2);

%% =========================================================
% 3. 扫频设置
%% =========================================================
Omega_Start = 10.0;
Omega_End   = 0.1;
Omega_Step  = -0.01;
Omega_Next  = Omega_Start + Omega_Step;

nStepsArc = 5000;

ParamMin = Omega_End;
ParamMax = Omega_Start + 0.05;

% 第二峰搜索频带
% 如果你的 baseline 第二峰不在这个范围内，改这里即可。
second_peak_band = [1.2, 4.0];

%% =========================================================
% 4. 先计算 baseline，用于自动识别第二峰
%% =========================================================
fprintf('\n====================================================\n');
fprintf('Step 1: compute baseline, no circuit.\n');

case_baseline.lam     = 0.0;
case_baseline.kap_e   = 0.0;
case_baseline.kap_c   = 0.0;
case_baseline.sigma   = 0.0;
case_baseline.Omega_e = NaN;
case_baseline.Csh     = NaN;
case_baseline.Lsh     = NaN;
case_baseline.Rt      = NaN;

sysP_baseline = make_sysP(alpha1, alpha2, gamma1, gamma2, mu, ...
                          case_baseline.lam, ...
                          case_baseline.kap_e, ...
                          case_baseline.kap_c, ...
                          case_baseline.sigma);

FixedOmega = [];
ParamMin = Omega_End;
ParamMax = Omega_Start + 0.05;

[x_base, info_base] = run_backward_frf(sysP_baseline, ...
                                       Omega_Start, ...
                                       Omega_Next, ...
                                       nStepsArc);

[Om_base, TF_base_dB, TF_base_lin] = calc_TF_from_branch(x_base, sysP_baseline, Fw);

valid = isfinite(Om_base) & isfinite(TF_base_dB) & Om_base > 0;
Om_base = Om_base(valid);
TF_base_dB = TF_base_dB(valid);
TF_base_lin = TF_base_lin(valid);

% 按 Omega 升序排序，方便找峰
[Om_base_sort, idx_sort] = sort(Om_base);
TF_base_dB_sort = TF_base_dB(idx_sort);

idx_second = Om_base_sort >= second_peak_band(1) & Om_base_sort <= second_peak_band(2);

if ~any(idx_second)
    error('No baseline points in second_peak_band = [%.3f, %.3f].', ...
          second_peak_band(1), second_peak_band(2));
end

Om_second_band = Om_base_sort(idx_second);
TF_second_band = TF_base_dB_sort(idx_second);

% 找该频带内最大值，作为第二峰
[baseline_peak2_dB, idx_p2] = max(TF_second_band);
baseline_peak2_Om = Om_second_band(idx_p2);

fprintf('\n===== Baseline second peak detected =====\n');
fprintf('Second peak search band: [%.3f, %.3f]\n', ...
        second_peak_band(1), second_peak_band(2));
fprintf('Baseline second peak: Omega = %.6f, TF = %.6f dB\n', ...
        baseline_peak2_Om, baseline_peak2_dB);

%% =========================================================
% 5. 根据 baseline 第二峰自动生成目标 Omega_e
%% =========================================================
% 目标：一个偏低、一个对准、两个偏高。
% 这样可以验证：Omega_e 接近第二峰时是否压制最好。
target_ratio = [0.75, 1.00, 1.25, 1.50];

Omega_e_list = baseline_peak2_Om * target_ratio;

% 也可以手动指定，例如：
% Omega_e_list = [1.6, 2.0, 2.5, 3.0];

fprintf('\n===== Target Omega_e for second peak test =====\n');
for i = 1:numel(Omega_e_list)
    fprintf('Omega_e_%d = %.6f\n', i, Omega_e_list(i));
end

%% =========================================================
% 6. 构造只改变 Csh 的电路工况
%% =========================================================
case_names = cell(1, numel(Omega_e_list)+1);
case_data  = struct([]);

case_names{1} = 'Baseline: no circuit';

case_data(1).lam       = 0.0;
case_data(1).kap_e     = 0.0;
case_data(1).kap_c     = 0.0;
case_data(1).sigma     = 0.0;
case_data(1).Omega_e   = NaN;
case_data(1).Csh       = NaN;
case_data(1).Lsh       = NaN;
case_data(1).Rt        = NaN;

for i = 1:numel(Omega_e_list)

    Omega_e_target = Omega_e_list(i);

    % 只改变电容 Csh，使 Omega_e 达到目标值
    Csh_i = 1/(Lsh * wn^2 * Omega_e_target^2);

    % 无量纲倒电容
    kap_c_i = 1/(Csh_i * R0 * wn);

    % 检查实际得到的无量纲电路固有频率
    Omega_e_check = sqrt(kap_c_i / kap_e_phys);

    case_names{i+1} = sprintf('\\Omega_e=%.3f, C=%.5f F', ...
                               Omega_e_check, Csh_i);

    case_data(i+1).lam       = lambda_phys;
    case_data(i+1).kap_e     = kap_e_phys;
    case_data(i+1).kap_c     = kap_c_i;
    case_data(i+1).sigma     = sigma_phys;
    case_data(i+1).Omega_e   = Omega_e_check;
    case_data(i+1).Csh       = Csh_i;
    case_data(i+1).Lsh       = Lsh;
    case_data(i+1).Rt        = Rt;
end
%% =========================================================
% 7. 循环计算所有工况
%% =========================================================
results = struct([]);

% 先保存 baseline
results(1).name = case_names{1};
results(1).failed = false;
results(1).case = case_data(1);
results(1).sysP = sysP_baseline;
results(1).x_res = x_base;
results(1).Om = Om_base;
results(1).TF_dB = TF_base_dB;
results(1).TF_lin = TF_base_lin;
results(1).info = info_base;

% baseline 第二峰统计
[results(1).peak2_dB, results(1).peak2_Om] = ...
    get_peak_in_band(Om_base, TF_base_dB, second_peak_band);

fprintf('\n====================================================\n');
fprintf('Baseline saved. Start Csh-tuned cases.\n');

for ic = 2:numel(case_data)

    fprintf('\n====================================================\n');
    fprintf('Case %d/%d: %s\n', ic, numel(case_data), case_names{ic});
    fprintf('lambda = %.8f, kap_e = %.8f, kap_c = %.8f, sigma = %.8f\n', ...
        case_data(ic).lam, ...
        case_data(ic).kap_e, ...
        case_data(ic).kap_c, ...
        case_data(ic).sigma);

    fprintf('Physical Csh = %.8f F, target/check Omega_e = %.6f\n', ...
        case_data(ic).Csh, case_data(ic).Omega_e);

    FixedOmega = [];
    ParamMin = Omega_End;
    ParamMax = Omega_Start + 0.05;

    sysP = make_sysP(alpha1, alpha2, gamma1, gamma2, mu, ...
                     case_data(ic).lam, ...
                     case_data(ic).kap_e, ...
                     case_data(ic).kap_c, ...
                     case_data(ic).sigma);

    try
        [x_res, info] = run_backward_frf(sysP, ...
                                         Omega_Start, ...
                                         Omega_Next, ...
                                         nStepsArc);
    catch ME
        warning('Case %d failed: %s', ic, ME.message);

        results(ic).name = case_names{ic};
        results(ic).failed = true;
        results(ic).message = ME.message;
        continue;
    end

    [Om, TF_dB, TF_lin] = calc_TF_from_branch(x_res, sysP, Fw);

    valid = isfinite(Om) & isfinite(TF_dB) & Om > 0;
    Om = Om(valid);
    TF_dB = TF_dB(valid);
    TF_lin = TF_lin(valid);

    results(ic).name = case_names{ic};
    results(ic).failed = false;
    results(ic).case = case_data(ic);
    results(ic).sysP = sysP;
    results(ic).x_res = x_res;
    results(ic).Om = Om;
    results(ic).TF_dB = TF_dB;
    results(ic).TF_lin = TF_lin;
    results(ic).info = info;

    [peak_val, peak_om] = get_peak_in_band(Om, TF_dB, second_peak_band);

    results(ic).peak2_dB = peak_val;
    results(ic).peak2_Om = peak_om;

    fprintf('Peak in [%.2f, %.2f]: TFmax = %.4f dB at Omega = %.5f\n', ...
        second_peak_band(1), second_peak_band(2), peak_val, peak_om);
end

%% =========================================================
% 8. FRF 全频段对比图
%% =========================================================
figure('Color','w', 'Position',[100 100 900 580]);
ax = gca; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
set(ax,'XScale','log');

for ic = 1:numel(results)
    if isfield(results(ic),'failed') && results(ic).failed
        continue;
    end

    if ic == 1
        plot(ax, results(ic).Om, results(ic).TF_dB, ...
            'k-', 'LineWidth', 2.2, ...
            'DisplayName', results(ic).name);
    else
        plot(ax, results(ic).Om, results(ic).TF_dB, ...
            'LineWidth', 1.6, ...
            'DisplayName', results(ic).name);
    end
end

yline(ax, 0, 'k--', '0 dB', 'HandleVisibility','off');

xlabel(ax, '\Omega');
ylabel(ax, 'Force Transmissibility 20log_{10}(|f_t|/f) (dB)');
title(ax, 'Second-peak tuning by changing C_{sh} only');

xlim(ax, [Omega_End, Omega_Start]);
ylim(ax, [-60, 20]);
legend(ax, 'Location','best');

% 标出目标频率和 baseline 第二峰
xline(ax, baseline_peak2_Om, 'k--', ...
    sprintf('Baseline peak2=%.3f', baseline_peak2_Om), ...
    'HandleVisibility','off');

for i = 1:numel(Omega_e_list)
    xline(ax, Omega_e_list(i), ':', ...
        sprintf('\\Omega_e=%.3f', Omega_e_list(i)), ...
        'HandleVisibility','off');
end

%% =========================================================
% 9. 第二峰区域局部放大图
%% =========================================================
figure('Color','w', 'Position',[150 150 900 580]);
ax2 = gca; hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
set(ax2,'XScale','log');

for ic = 1:numel(results)
    if isfield(results(ic),'failed') && results(ic).failed
        continue;
    end

    if ic == 1
        plot(ax2, results(ic).Om, results(ic).TF_dB, ...
            'k-', 'LineWidth', 2.2, ...
            'DisplayName', results(ic).name);
    else
        plot(ax2, results(ic).Om, results(ic).TF_dB, ...
            'LineWidth', 1.6, ...
            'DisplayName', results(ic).name);
    end
end

yline(ax2, 0, 'k--', '0 dB', 'HandleVisibility','off');

xlabel(ax2, '\Omega');
ylabel(ax2, 'Force Transmissibility 20log_{10}(|f_t|/f) (dB)');
title(ax2, 'Zoomed view near the second peak');

xlim(ax2, second_peak_band);

% 自动设置 y 轴范围
all_y = [];
for ic = 1:numel(results)
    if isfield(results(ic),'failed') && results(ic).failed
        continue;
    end
    idx_tmp = results(ic).Om >= second_peak_band(1) & ...
              results(ic).Om <= second_peak_band(2);
    all_y = [all_y; results(ic).TF_dB(idx_tmp)];
end

if ~isempty(all_y)
    ylim(ax2, [min(all_y)-5, max(all_y)+5]);
end

legend(ax2, 'Location','best');

xline(ax2, baseline_peak2_Om, 'k--', ...
    sprintf('Baseline peak2=%.3f', baseline_peak2_Om), ...
    'HandleVisibility','off');

for i = 1:numel(Omega_e_list)
    xline(ax2, Omega_e_list(i), ':', ...
        sprintf('\\Omega_e=%.3f', Omega_e_list(i)), ...
        'HandleVisibility','off');
end

%% =========================================================
% 10. Baseline 第二峰与归一化 c_eq 对齐图
%% =========================================================
Omega_grid = logspace(log10(Omega_End), log10(Omega_Start), 1500);

figure('Color','w', 'Position',[200 200 900 580]);

yyaxis left
hold on; grid on; box on;
set(gca,'XScale','log');

plot(Om_base, TF_base_dB, 'k-', 'LineWidth', 2.2, ...
    'DisplayName', 'Baseline TF');

ylabel('Baseline TF (dB)');
ylim([-60, 20]);
yline(0, 'k--', '0 dB', 'HandleVisibility','off');

yyaxis right
hold on;

for ic = 2:numel(case_data)

    lambda = case_data(ic).lam;
    kap_e  = case_data(ic).kap_e;
    kap_c  = case_data(ic).kap_c;
    sigma  = case_data(ic).sigma;

    [~, ~, ceq] = complex_operator_parts(Omega_grid, ...
                                         lambda, ...
                                         kap_e, ...
                                         kap_c, ...
                                         sigma);

    ceq_norm = ceq ./ max(abs(ceq));

    semilogx(Omega_grid, ceq_norm, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('c_{eq}, \\Omega_e=%.3f', ...
        case_data(ic).Omega_e));
end

ylabel('Normalized c_{eq}(\Omega)');
ylim([0, 1.1]);

xlabel('\Omega');
title('Baseline second peak and equivalent damping-shaping regions');
xlim([Omega_End, Omega_Start]);
legend('Location','best');

xline(baseline_peak2_Om, 'k--', ...
      sprintf('Baseline peak2=%.3f', baseline_peak2_Om), ...
      'HandleVisibility','off');

%% =========================================================
% 11. 输出第二峰压制效果表
%% =========================================================
fprintf('\n\n====================================================\n');
fprintf('Second-peak comparison in Omega band [%.2f, %.2f]\n', ...
        second_peak_band(1), second_peak_band(2));
fprintf('Only Csh is changed. Kt, Ke, R0, Rt, Lsh, m1, k1 fixed.\n');
fprintf('====================================================\n');

baseline_peak = results(1).peak2_dB;

fprintf('%-22s %-14s %-14s %-14s %-14s %-14s\n', ...
    'Case', 'Csh(F)', 'Omega_peak', 'Peak_dB', 'Reduction_dB', 'Omega_e');

for ic = 1:numel(results)

    if isfield(results(ic),'failed') && results(ic).failed
        fprintf('%-22s %-14s %-14s %-14s %-14s %-14s\n', ...
            results(ic).name, 'FAILED', '-', '-', '-', '-');
        continue;
    end

    if ic == 1
        fprintf('%-22s %-14s %-14.5f %-14.4f %-14s %-14s\n', ...
            'Baseline', ...
            '-', ...
            results(ic).peak2_Om, ...
            results(ic).peak2_dB, ...
            '-', ...
            '-');
    else
        reduction = baseline_peak - results(ic).peak2_dB;

        fprintf('%-22s %-14.8f %-14.5f %-14.4f %-14.4f %-14.5f\n', ...
            sprintf('C-tuned %d', ic-1), ...
            results(ic).case.Csh, ...
            results(ic).peak2_Om, ...
            results(ic).peak2_dB, ...
            reduction, ...
            results(ic).case.Omega_e);
    end
end

fprintf('\n判断逻辑：\n');
fprintf('如果 Omega_e 接近 baseline 第二峰频率时，第二峰频带内 Peak_dB 降低最多，\n');
fprintf('就说明复算子阻尼整形频带对准第二峰的设计逻辑成立。\n');
fprintf('如果效果不明显，说明该峰可能主要受高频附加刚度、模态耦合或非线性分支影响，\n');
fprintf('单纯调 Csh 对第二峰的控制能力有限，需要进一步调 Rt 或 Lsh。\n');

%% =========================================================
% 12. 诊断图：第二峰附近的原始 Kr 和原始 ceq
%     同色对应：
%       实线  = K_r
%       虚线  = c_eq
%       点线竖线 = 对应 Omega_e
%% =========================================================

Omega_grid2 = linspace(second_peak_band(1), second_peak_band(2), 1000);

% 有效电路工况数量，不含 baseline
nCircuitCases = numel(case_data) - 1;

% 为每个 Omega_e 分配一种颜色
clr = lines(nCircuitCases);

figure('Color','w', 'Position',[250 250 950 600]);
ax = gca;
hold(ax, 'on');
grid(ax, 'on');
box(ax, 'on');

%% ---------- 左轴：K_r ----------
yyaxis left
hold on;

hKr = gobjects(nCircuitCases, 1);

for ii = 1:nCircuitCases

    ic = ii + 1;   % case_data(1) 是 baseline，所以电路工况从 2 开始

    lambda = case_data(ic).lam;
    kap_e  = case_data(ic).kap_e;
    kap_c  = case_data(ic).kap_c;
    sigma  = case_data(ic).sigma;

    [Kr, ~, ~] = complex_operator_parts(Omega_grid2, ...
                                         lambda, ...
                                         kap_e, ...
                                         kap_c, ...
                                         sigma);

    hKr(ii) = plot(Omega_grid2, Kr, '-', ...
        'Color', clr(ii,:), ...
        'LineWidth', 2.0, ...
        'DisplayName', sprintf('K_r, \\Omega_e=%.3f', case_data(ic).Omega_e));
end

ylabel('K_r(\Omega) = k_{eq}(\Omega)');
yline(0, 'k:', 'LineWidth', 1.0, 'HandleVisibility','off');

%% ---------- 右轴：c_eq ----------
yyaxis right
hold on;

hCeq = gobjects(nCircuitCases, 1);

for ii = 1:nCircuitCases

    ic = ii + 1;

    lambda = case_data(ic).lam;
    kap_e  = case_data(ic).kap_e;
    kap_c  = case_data(ic).kap_c;
    sigma  = case_data(ic).sigma;

    [~, ~, ceq] = complex_operator_parts(Omega_grid2, ...
                                         lambda, ...
                                         kap_e, ...
                                         kap_c, ...
                                         sigma);

    hCeq(ii) = plot(Omega_grid2, ceq, '--', ...
        'Color', clr(ii,:), ...
        'LineWidth', 2.0, ...
        'DisplayName', sprintf('c_{eq}, \\Omega_e=%.3f', case_data(ic).Omega_e));
end

ylabel('c_{eq}(\Omega), unnormalized');

%% ---------- 公共设置 ----------
xlabel('\Omega');
title('Second-peak diagnostic: K_r and c_{eq} with matched colors');
xlim(second_peak_band);

% baseline 第二峰位置
xline(baseline_peak2_Om, 'k--', ...
      sprintf('Baseline peak2=%.3f', baseline_peak2_Om), ...
      'LineWidth', 1.2, ...
      'HandleVisibility','off');

% 对应每组电气固有频率 Omega_e，用同色点线标出
for ii = 1:nCircuitCases

    ic = ii + 1;

    xline(case_data(ic).Omega_e, ':', ...
        sprintf('\\Omega_e=%.3f', case_data(ic).Omega_e), ...
        'Color', clr(ii,:), ...
        'LineWidth', 1.3, ...
        'HandleVisibility','off');
end

%% ---------- 图例 ----------
% 图例中保留每条曲线：
% 同色实线/虚线分别对应同一 Omega_e 下的 Kr 和 ceq
legend([hKr; hCeq], 'Location','best');

%% ---------- 添加说明文字 ----------
annotation('textbox', [0.15 0.82 0.35 0.08], ...
    'String', {'Same color = same \Omega_e / same C_{sh}', ...
               'Solid line: K_r(\Omega), dashed line: c_{eq}(\Omega)'}, ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'w', ...
    'EdgeColor', [0.7 0.7 0.7]);

%% =========================================================
% 13. 输出第二峰位置处的 Kr 和 ceq
%% =========================================================
fprintf('\n\n====================================================\n');
fprintf('Operator values at baseline second peak Omega = %.6f\n', baseline_peak2_Om);
fprintf('====================================================\n');

fprintf('%-14s %-14s %-14s %-14s %-14s\n', ...
    'Omega_e', 'Csh(F)', 'Kr@peak2', 'ceq@peak2', 'Peak_dB');

for ic = 2:numel(case_data)

    lambda = case_data(ic).lam;
    kap_e  = case_data(ic).kap_e;
    kap_c  = case_data(ic).kap_c;
    sigma  = case_data(ic).sigma;

    [Kr_p, ~, ceq_p] = complex_operator_parts(baseline_peak2_Om, ...
                                              lambda, ...
                                              kap_e, ...
                                              kap_c, ...
                                              sigma);

    fprintf('%-14.6f %-14.8f %-14.6f %-14.6f %-14.4f\n', ...
        case_data(ic).Omega_e, ...
        case_data(ic).Csh, ...
        Kr_p, ...
        ceq_p, ...
        results(ic).peak2_dB);
end
%% =========================================================
% Local functions
%% =========================================================

function sysP = make_sysP(alpha1, alpha2, gamma1, gamma2, mu, lam, kap_e, kap_c, sigma)

    % 你的代码中上层线性项采用 be1 + al1
    P.be1 = 1.0;
    P.al1 = alpha1 - P.be1;

    P.be2 = alpha2;
    P.ga1 = gamma1;
    P.ga2 = gamma2;

    P.mu  = mu;
    P.ze1 = 0.05;    % 下层对地阻尼 zeta2

    P.lam   = lam;
    P.kap_e = kap_e;
    P.kap_c = kap_c;
    P.sigma = sigma;

    sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, ...
            P.lam, P.kap_e, P.kap_c, P.sigma, P.ga2];
end


function [x_res, info] = run_backward_frf(sysP, Omega_Start, Omega_Next, nStepsArc)

    global FixedOmega ParamMin ParamMax

    FixedOmega = [];

    if isempty(ParamMin)
        ParamMin = 0.1;
    end

    if isempty(ParamMax)
        ParamMax = Omega_Start + 0.05;
    end

    % 高频起点
    y_init = zeros(16,1);
    y_init(16) = Omega_Start;

    [x0_full, ok0, R0] = newton('nondim_temp2', y_init, sysP);

    if ~ok0 || R0 > 1e-6
        error('High-frequency starting point failed: Omega=%.6f, R=%.3e', ...
              Omega_Start, R0);
    end

    x0 = x0_full(1:15);

    % 第二个点
    y_init2 = [x0; Omega_Next];

    [x1_full, ok1, R1] = newton('nondim_temp2', y_init2, sysP);

    if ~ok1 || R1 > 1e-6
        error('Second point failed: Omega=%.6f, R=%.3e', ...
              Omega_Next, R1);
    end

    x1 = x1_full(1:15);

    % 弧长延拓
    [x_res, info] = branch_follow2('nondim_temp2', ...
                                   nStepsArc, ...
                                   Omega_Start, ...
                                   Omega_Next, ...
                                   x0, x1, ...
                                   sysP);
end


function [Om, TF_dB, TF_lin] = calc_TF_from_branch(x_res, sysP, Fw)

    Om  = x_res(16,:).';

    be2 = sysP(2);
    mu  = sysP(3);
    ze2 = sysP(6);
    ga2 = sysP(11);

    % 下层位移 x2 的 0/1/3 谐波系数
    x2 = x_res(6:10,:).';

    % x2_dot
    W = Om;
    x2_dot = zeros(size(x2));

    x2_dot(:,1) = 0;
    x2_dot(:,2) = W .* x2(:,3);
    x2_dot(:,3) = -W .* x2(:,2);
    x2_dot(:,4) = 3*W .* x2(:,5);
    x2_dot(:,5) = -3*W .* x2(:,4);

    % x2^3 的 AFT 投影
    x2_cub = cubic_proj_013_batch_local(x2);

    % 基础传递力：
    % ft = alpha2*x2 + gamma2*x2^3 + 2*mu*zeta2*x2'
    ft = be2*x2 + ga2*x2_cub + 2*mu*ze2*x2_dot;

    ft1 = hypot(ft(:,2), ft(:,3));
    ft3 = hypot(ft(:,4), ft(:,5));

    % 基波与三次谐波合成幅值
    ft_amp = hypot(ft1, ft3);

    TF_lin = ft_amp ./ Fw;
    TF_dB  = 20*log10(max(TF_lin, 1e-300));
end


function cubic = cubic_proj_013_batch_local(U)

    [~, T_mat, T_inv] = get_AFT_matrices_local();

    X_time  = (T_mat * U.').';
    X3_time = X_time.^3;

    cubic = (T_inv * X3_time.').';
end


function [N, T_mat, T_inv] = get_AFT_matrices_local()

    persistent pN pT pTinv

    if isempty(pN)

        pN = 64;

        t = (0:pN-1)'*(2*pi/pN);

        c1 = cos(t);
        s1 = sin(t);
        c3 = cos(3*t);
        s3 = sin(3*t);
        dc = ones(pN,1);

        pT = [dc, c1, s1, c3, s3];

        Inv = [dc, 2*c1, 2*s1, 2*c3, 2*s3]';
        pTinv = (1/pN) * Inv;
        pTinv(1,:) = (1/pN) * dc';
    end

    N = pN;
    T_mat = pT;
    T_inv = pTinv;
end


function [Kr, Ki, ceq] = complex_operator_parts(Omega, lambda, kap_e, kap_c, sigma)

    % K = Kr + j Ki
    %
    % Kr = lambda*Omega^2*(kap_e*Omega^2-kap_c) / Den
    % Ki = lambda*sigma*Omega^3 / Den
    % ceq = Ki/Omega

    Den = (kap_e.*Omega.^2 - kap_c).^2 + (sigma.*Omega).^2;

    Kr = lambda .* Omega.^2 .* (kap_e.*Omega.^2 - kap_c) ./ Den;
    Ki = lambda .* sigma .* Omega.^3 ./ Den;

    ceq = Ki ./ Omega;
end


function [peak_val, peak_om] = get_peak_in_band(Om, TF_dB, band)

    valid = isfinite(Om) & isfinite(TF_dB) & ...
            Om >= band(1) & Om <= band(2);

    if ~any(valid)
        peak_val = NaN;
        peak_om = NaN;
        return;
    end

    Om_band = Om(valid);
    TF_band = TF_dB(valid);

    [Om_band, idx] = sort(Om_band);
    TF_band = TF_band(idx);

    [peak_val, idx_peak] = max(TF_band);
    peak_om = Om_band(idx_peak);
end