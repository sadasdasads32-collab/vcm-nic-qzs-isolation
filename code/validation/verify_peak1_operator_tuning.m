%% =========================================================
% verify_peak1_by_changing_Csh.m
%
% 目的：
%   验证第一个主峰附近，电路固有频率 Omega_e 的匹配效果。
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
%   只改变外接电容 Csh，使电路固有频率取：
%
%       Omega_e = 0.224, 0.300, 0.387, 0.500
%
%   其中：
%
%       Omega_e = 1/(wn*sqrt(Lsh*Csh))
%
%   所以：
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
% 3. 目标电路固有频率：只通过改变 Csh 实现
%% =========================================================
Omega_e_list = [0.224, 0.300, 0.387, 0.500];

% 第 1 个 case 是 baseline，不接电路
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

    case_names{i+1} = sprintf('\\Omega_e=%.3f, C=%.4f F', ...
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

fprintf('\n===== Target Omega_e by changing Csh only =====\n');
fprintf('%-12s %-14s %-14s %-14s %-14s\n', ...
        'Omega_e', 'Csh(F)', 'kap_e', 'kap_c', 'sigma');

for i = 2:numel(case_data)
    fprintf('%-12.6f %-14.8f %-14.8f %-14.8f %-14.8f\n', ...
        case_data(i).Omega_e, ...
        case_data(i).Csh, ...
        case_data(i).kap_e, ...
        case_data(i).kap_c, ...
        case_data(i).sigma);
end

%% =========================================================
% 4. 扫频设置
%% =========================================================
Omega_Start = 10.0;
Omega_End   = 0.1;
Omega_Step  = -0.01;
Omega_Next  = Omega_Start + Omega_Step;

nStepsArc = 5000;

ParamMin = Omega_End;
ParamMax = Omega_Start + 0.05;

%% =========================================================
% 5. 循环计算所有工况
%% =========================================================
results = struct([]);

for ic = 1:numel(case_data)

    fprintf('\n====================================================\n');
    fprintf('Case %d/%d: %s\n', ic, numel(case_data), case_names{ic});
    fprintf('lambda = %.8f, kap_e = %.8f, kap_c = %.8f, sigma = %.8f\n', ...
        case_data(ic).lam, ...
        case_data(ic).kap_e, ...
        case_data(ic).kap_c, ...
        case_data(ic).sigma);

    if ic > 1
        fprintf('Physical Csh = %.8f F, target/check Omega_e = %.6f\n', ...
            case_data(ic).Csh, case_data(ic).Omega_e);
    end

    % 每个 case 重新设置全局扫频状态，防止污染
    FixedOmega = [];
    ParamMin = Omega_End;
    ParamMax = Omega_Start + 0.05;

    % 组装 sysP
    sysP = make_sysP(alpha1, alpha2, gamma1, gamma2, mu, ...
                     case_data(ic).lam, ...
                     case_data(ic).kap_e, ...
                     case_data(ic).kap_c, ...
                     case_data(ic).sigma);

    % 跑 FRF
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

    % 计算 TF
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

    % 第一个主峰频段统计
    peak_band = [0.15, 0.60];
    idx_band = Om >= peak_band(1) & Om <= peak_band(2);

    if any(idx_band)
        Om_band = Om(idx_band);
        TF_band = TF_dB(idx_band);

        [peak_val, idx_local] = max(TF_band);
        peak_om = Om_band(idx_local);

        results(ic).peak1_dB = peak_val;
        results(ic).peak1_Om = peak_om;

        fprintf('Peak in [%.2f, %.2f]: TFmax = %.4f dB at Omega = %.5f\n', ...
            peak_band(1), peak_band(2), peak_val, peak_om);
    else
        results(ic).peak1_dB = NaN;
        results(ic).peak1_Om = NaN;

        fprintf('Warning: no valid point in first-peak band.\n');
    end
end

%% =========================================================
% 6. FRF 对比图
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
title(ax, 'First-peak tuning by changing C_{sh} only');

xlim(ax, [Omega_End, Omega_Start]);
ylim(ax, [-60, 20]);
legend(ax, 'Location','best');

% 标出目标频率
for i = 1:numel(Omega_e_list)
    xline(ax, Omega_e_list(i), ':', ...
        sprintf('\\Omega_e=%.3f', Omega_e_list(i)), ...
        'HandleVisibility','off');
end

%% =========================================================
% 7. Baseline TF 与等效阻尼 c_eq 对齐图
%% =========================================================
Omega_grid = logspace(log10(Omega_End), log10(Omega_Start), 1500);

figure('Color','w', 'Position',[150 150 900 580]);

yyaxis left
hold on; grid on; box on;
set(gca,'XScale','log');

if ~results(1).failed
    plot(results(1).Om, results(1).TF_dB, 'k-', 'LineWidth', 2.2, ...
        'DisplayName', 'Baseline TF');
end

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
title('Baseline first peak and equivalent damping-shaping regions');
xlim([Omega_End, Omega_Start]);
legend('Location','best');

%% =========================================================
% 8. 输出第一主峰压制效果表
%% =========================================================
fprintf('\n\n====================================================\n');
fprintf('First-peak comparison in Omega band [0.15, 0.60]\n');
fprintf('Only Csh is changed. Kt, Ke, R0, Rt, Lsh, m1, k1 fixed.\n');
fprintf('====================================================\n');

valid_baseline = isfield(results(1),'peak1_dB') && ~isnan(results(1).peak1_dB);
if valid_baseline
    baseline_peak = results(1).peak1_dB;
else
    baseline_peak = NaN;
end

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
            results(ic).peak1_Om, ...
            results(ic).peak1_dB, ...
            '-', ...
            '-');
    else
        reduction = baseline_peak - results(ic).peak1_dB;

        fprintf('%-22s %-14.8f %-14.5f %-14.4f %-14.4f %-14.5f\n', ...
            sprintf('C-tuned %d', ic-1), ...
            results(ic).case.Csh, ...
            results(ic).peak1_Om, ...
            results(ic).peak1_dB, ...
            reduction, ...
            results(ic).case.Omega_e);
    end
end

fprintf('\n判断逻辑：\n');
fprintf('如果 Omega_e 接近 baseline 第一主峰频率时，Peak_dB 降低最多，\n');
fprintf('就说明复算子阻尼整形频带对准主峰的设计逻辑成立。\n');


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

    if nargin('newton') < 3
        % 防御性分支，一般不会走到这里
    end

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