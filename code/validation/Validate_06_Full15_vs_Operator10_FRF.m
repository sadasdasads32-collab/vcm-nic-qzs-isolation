%% =========================================================
% Validate_06_Full15_vs_Operator10_FRF_FIXED.m
%
% 目的：
% 验证 15 维完整机电模型与 10 维复算子消元模型的 FRF 等价性。
%
% 说明：
% 1) Part D-1 保持延拓顺序画图，用于观察多值分支；
% 2) Part D-2 排序 + 插值，只用于误差统计；
% 3) 兼容新版 branch_follow2 / branch_follow2N 返回 info 结构体。
% =========================================================

clc; clear; close all;

%% ---------------------------------------------------------
% 0) 全局变量
%% ---------------------------------------------------------
global Fw FixedOmega ParamMin ParamMax

Fw = 0.005;
FixedOmega = [];

ParamMin = 0.1;
ParamMax = 10.0;

%% ---------------------------------------------------------
% 1) 物理电路参数
%% ---------------------------------------------------------
Kt = 7.474;
Ke = 7.474;
m1 = 2.2;
k1 = 3000;
R0 = 3.8;

wn = sqrt(k1/m1);

Rt  = 2.3674;
Lsh = 0.04065;
Csh = 0.2227;

P.lam   = Kt*Ke*wn/(k1*R0);
P.kap_e = Lsh*wn/R0;
P.kap_c = 1/(Csh*R0*wn);
P.sigma = Rt/R0;

%% ---------------------------------------------------------
% 2) 机械/QZS 参数
%% ---------------------------------------------------------
mu   = 0.2;
beta = 2.0;
K1   = 1.0;
K2   = 0.0;
U    = 2.0;
L    = 4/9;

v = 2.5;

alpha1 = v    - 2*K1*(1-L)/L;
alpha2 = beta - 2*K2*(1-L)/L;
gamma1 = K1/(U^2 * L^3);
gamma2 = K2/(U^2 * L^3);

P.be1 = 1.0;
P.al1 = alpha1 - P.be1;
P.be2 = alpha2;
P.ga1 = gamma1;
P.ga2 = gamma2;
P.mu  = mu;
P.ze1 = 0.05;

sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, ...
        P.lam, P.kap_e, P.kap_c, P.sigma, P.ga2];

fprintf('\n========== 参数检查 ==========\n');
fprintf('lambda  = %.8f\n', P.lam);
fprintf('theta   = %.8f\n', sqrt(P.lam));
fprintf('kap_e   = %.8f\n', P.kap_e);
fprintf('kap_c   = %.8f\n', P.kap_c);
fprintf('sigma   = %.8f\n', P.sigma);
fprintf('alpha1  = %.8f\n', P.be1 + P.al1);
fprintf('alpha2  = %.8f\n', P.be2);
fprintf('gamma1  = %.8f\n', P.ga1);
fprintf('gamma2  = %.8f\n', P.ga2);

%% ---------------------------------------------------------
% 3) 延拓设置
%% ---------------------------------------------------------
Omega_Start = 10.0;
Omega_Step  = -0.01;
Omega_Next  = Omega_Start + Omega_Step;

Nsteps_full = 3000;
Nsteps_op   = 3000;

fprintf('\n========== 延拓设置 ==========\n');
fprintf('Omega_Start = %.6f\n', Omega_Start);
fprintf('Omega_Next  = %.6f\n', Omega_Next);
fprintf('ParamMin    = %.6f\n', ParamMin);
fprintf('ParamMax    = %.6f\n', ParamMax);

%% =========================================================
% Part A: 15维完整机电模型
%% =========================================================
fprintf('\n\n=================================================\n');
fprintf('开始计算 15 维完整机电模型 FRF\n');
fprintf('=================================================\n');

%% A1) 初始两点
y0_full = zeros(16,1);
y0_full(16) = Omega_Start;

[x0_full_all, ok0_full, R0_full] = newton('nondim_temp2', y0_full, sysP);

fprintf('15D point0: Omega=%.6f, ok=%d, R=%.3e\n', ...
        Omega_Start, ok0_full, R0_full);

if ~ok0_full || R0_full > 1e-6
    error('15D: 高频起点 Newton 失败或残差过大。');
end

x0_full = x0_full_all(1:15);

y1_full = [x0_full; Omega_Next];

[x1_full_all, ok1_full, R1_full] = newton('nondim_temp2', y1_full, sysP);

fprintf('15D point1: Omega=%.6f, ok=%d, R=%.3e\n', ...
        Omega_Next, ok1_full, R1_full);

if ~ok1_full || R1_full > 1e-6
    error('15D: 第二点 Newton 失败或残差过大。');
end

x1_full = x1_full_all(1:15);

%% A2) 弧长延拓
[xres_full, info_full] = branch_follow2('nondim_temp2', ...
                                        Nsteps_full, ...
                                        Omega_Start, Omega_Next, ...
                                        x0_full, x1_full, sysP);

msg_full = get_stop_message(info_full);

fprintf('15D continuation done. stop=%s, points=%d\n', ...
        msg_full, size(xres_full,2));

%% A3) 后处理 TF
Om_full_raw = xres_full(16,:).';
TF_full_raw = post_TF_full(xres_full, sysP, Fw);

valid_full = isfinite(Om_full_raw) & isfinite(TF_full_raw) & Om_full_raw > 0;

% 注意：这里保持延拓顺序，后面用来画多值分支
Om_full = Om_full_raw(valid_full);
TF_full_dB = TF_full_raw(valid_full);

%% =========================================================
% Part B: 10维复算子消元模型
%% =========================================================
fprintf('\n\n=================================================\n');
fprintf('开始计算 10 维复算子消元模型 FRF\n');
fprintf('=================================================\n');

%% B1) 初始两点
y0_op = zeros(11,1);
y0_op(11) = Omega_Start;

[x0_op_all, ok0_op, R0_op] = newton('nondim_temp2_op', y0_op, sysP);

fprintf('10D point0: Omega=%.6f, ok=%d, R=%.3e\n', ...
        Omega_Start, ok0_op, R0_op);

if ~ok0_op || R0_op > 1e-6
    error('10D: 高频起点 Newton 失败或残差过大。');
end

x0_op = x0_op_all(1:10);

y1_op = [x0_op; Omega_Next];

[x1_op_all, ok1_op, R1_op] = newton('nondim_temp2_op', y1_op, sysP);

fprintf('10D point1: Omega=%.6f, ok=%d, R=%.3e\n', ...
        Omega_Next, ok1_op, R1_op);

if ~ok1_op || R1_op > 1e-6
    error('10D: 第二点 Newton 失败或残差过大。');
end

x1_op = x1_op_all(1:10);

%% B2) 弧长延拓
[xres_op, info_op] = branch_follow2N('nondim_temp2_op', ...
                                     Nsteps_op, ...
                                     Omega_Start, Omega_Next, ...
                                     x0_op, x1_op, sysP);

msg_op = get_stop_message(info_op);

fprintf('10D continuation done. stop=%s, points=%d\n', ...
        msg_op, size(xres_op,2));

%% B3) 后处理 TF
Om_op_raw = xres_op(11,:).';
TF_op_raw = post_TF_op(xres_op, sysP, Fw);

valid_op = isfinite(Om_op_raw) & isfinite(TF_op_raw) & Om_op_raw > 0;

% 注意：这里保持延拓顺序，后面用来画多值分支
Om_op = Om_op_raw(valid_op);
TF_op_dB = TF_op_raw(valid_op);

%% =========================================================
% Part C: 15D vs 10D 误差比较
% 注意：这里会排序和去重，只用于插值误差统计，不用于展示多值分支。
%% =========================================================
fprintf('\n\n=================================================\n');
fprintf('开始比较 15D 完整模型 vs 10D 算子模型\n');
fprintf('=================================================\n');

%% C1) 排序 + 去重
[Om_full_s, idxF] = sort(Om_full(:), 'ascend');
TF_full_s = TF_full_dB(idxF);

[Om_op_s, idxO] = sort(Om_op(:), 'ascend');
TF_op_s = TF_op_dB(idxO);

% 去除重复 Omega，避免 interp1 报错
[Om_full_u, iaF] = unique(Om_full_s, 'stable');
TF_full_u = TF_full_s(iaF);

[Om_op_u, iaO] = unique(Om_op_s, 'stable');
TF_op_u = TF_op_s(iaO);

%% C2) 共同频率范围
Om_min_cmp = max(min(Om_full_u), min(Om_op_u));
Om_max_cmp = min(max(Om_full_u), max(Om_op_u));

valid_cmp = Om_op_u >= Om_min_cmp & Om_op_u <= Om_max_cmp;

Om_cmp = Om_op_u(valid_cmp);
TF_op_cmp_dB = TF_op_u(valid_cmp);

% 将 15D 完整模型插值到 10D 的频率点
TF_full_cmp_dB = interp1(Om_full_u, TF_full_u, Om_cmp, 'linear');

diff_dB = TF_op_cmp_dB - TF_full_cmp_dB;

Einf_dB = max(abs(diff_dB));
Erms_dB = sqrt(mean(diff_dB.^2));

% 线性幅值域误差
TF_op_lin   = 10.^(TF_op_cmp_dB/20);
TF_full_lin = 10.^(TF_full_cmp_dB/20);

E2_lin = norm(TF_op_lin - TF_full_lin) / max(1e-14, norm(TF_full_lin));
Einf_lin = max(abs(TF_op_lin - TF_full_lin));

fprintf('\n========== 等价性误差指标 ==========\n');
fprintf('Compared points      = %d\n', numel(Om_cmp));
fprintf('Omega common range   = [%.6f, %.6f]\n', Om_min_cmp, Om_max_cmp);
fprintf('Max abs error dB     = %.6e dB\n', Einf_dB);
fprintf('RMS error dB         = %.6e dB\n', Erms_dB);
fprintf('Relative L2 error    = %.6e\n', E2_lin);
fprintf('Max linear TF error  = %.6e\n', Einf_lin);

% 更合理的工程判据
if Einf_dB < 5e-2 && Erms_dB < 1e-2 && E2_lin < 1e-3
    fprintf('判定：15D 完整模型与 10D 算子模型高度一致。\n');
elseif Einf_dB < 1e-1 && Erms_dB < 5e-2 && E2_lin < 5e-3
    fprintf('判定：15D 与 10D 基本一致，误差可接受。\n');
else
    fprintf('警告：15D 与 10D 差异偏大，请检查算子模型符号或后处理方式。\n');
end

%% =========================================================
% Part D-1: 按延拓顺序画图，用于观察多值分支
% 不排序，不 unique，不插值。
%% =========================================================
figure('Color','w','Position',[100 100 900 600]);
hold on; grid on; box on;

plot(Om_full, TF_full_dB, 'k.-', 'LineWidth', 1.2, 'MarkerSize', 7);
plot(Om_op,   TF_op_dB,   'r.--', 'LineWidth', 1.2, 'MarkerSize', 7);

xlabel('\Omega');
ylabel('TF (dB)');
title('15D full vs 10D operator: continuation-order plot');
legend('15D full','10D operator','Location','best');

xlim([ParamMin, ParamMax]);

%% =========================================================
% Part D-2: 排序后的等价性误差图
% 只用于单值插值比较，不用于展示多值分支。
%% =========================================================
figure('Color','w','Position',[100 100 1200 520]);

subplot(1,2,1);
hold on; grid on; box on;
set(gca,'XScale','log');

plot(Om_full_u, TF_full_u, 'k-', 'LineWidth', 2.0);
plot(Om_op_u,   TF_op_u,   'r--', 'LineWidth', 1.5);

xlabel('\Omega (log)');
ylabel('TF (dB)');
title('15D full model vs 10D operator model');
legend('15D full','10D operator','Location','best');
xlim([ParamMin, ParamMax]);

subplot(1,2,2);
hold on; grid on; box on;
set(gca,'XScale','log');

plot(Om_cmp, diff_dB, 'b-', 'LineWidth', 1.4);
yline(0, 'k--');

xlabel('\Omega (log)');
ylabel('\Delta TF = TF_{op} - TF_{full} (dB)');
title(sprintf('Difference, max=%.2e dB, RMS=%.2e dB', Einf_dB, Erms_dB));
xlim([Om_min_cmp, Om_max_cmp]);

%% 导出到工作区
assignin('base', 'Om_full_order', Om_full);
assignin('base', 'TF_full_dB_order', TF_full_dB);
assignin('base', 'Om_op_order', Om_op);
assignin('base', 'TF_op_dB_order', TF_op_dB);

assignin('base', 'Om_full_sorted', Om_full_u);
assignin('base', 'TF_full_dB_sorted', TF_full_u);
assignin('base', 'Om_op_sorted', Om_op_u);
assignin('base', 'TF_op_dB_sorted', TF_op_u);
assignin('base', 'Om_cmp', Om_cmp);
assignin('base', 'diff_dB', diff_dB);

fprintf('\n验证完成。\n');

%% =========================================================
% Helper 0: 兼容新版/旧版 continuation 输出信息
%% =========================================================
function msg = get_stop_message(info)
    if isstruct(info)
        if isfield(info, 'stop_reason')
            msg = info.stop_reason;
        else
            msg = 'info_struct';
        end
    else
        try
            msg = char(string(info));
        catch
            msg = 'unknown';
        end
    end
end

%% =========================================================
% Helper 1: 15D 完整模型 TF 后处理
%% =========================================================
function TF_dB = post_TF_full(xres, sysP, Fw)

    Om  = xres(16,:).';
    be2 = sysP(2);
    mu  = sysP(3);
    ze2 = sysP(6);
    ga2 = sysP(11);

    x2 = xres(6:10,:).';

    x2_dot = zeros(size(x2));
    x2_dot(:,1) = 0;
    x2_dot(:,2) = Om .* x2(:,3);
    x2_dot(:,3) = -Om .* x2(:,2);
    x2_dot(:,4) = 3*Om .* x2(:,5);
    x2_dot(:,5) = -3*Om .* x2(:,4);

    x2_cub = cubic_proj_013_batch(x2);

    ft = be2*x2 + ga2*x2_cub + 2*mu*ze2*x2_dot;

    ft1 = hypot(ft(:,2), ft(:,3));
    ft3 = hypot(ft(:,4), ft(:,5));
    ft_amp = hypot(ft1, ft3);

    TF_dB = 20*log10(max(ft_amp ./ Fw, 1e-300));
end

%% =========================================================
% Helper 2: 10D 算子模型 TF 后处理
%% =========================================================
function TF_dB = post_TF_op(xres, sysP, Fw)

    Om  = xres(11,:).';
    be2 = sysP(2);
    mu  = sysP(3);
    ze2 = sysP(6);
    ga2 = sysP(11);

    x2 = xres(6:10,:).';

    x2_dot = zeros(size(x2));
    x2_dot(:,1) = 0;
    x2_dot(:,2) = Om .* x2(:,3);
    x2_dot(:,3) = -Om .* x2(:,2);
    x2_dot(:,4) = 3*Om .* x2(:,5);
    x2_dot(:,5) = -3*Om .* x2(:,4);

    x2_cub = cubic_proj_013_batch(x2);

    ft = be2*x2 + ga2*x2_cub + 2*mu*ze2*x2_dot;

    ft1 = hypot(ft(:,2), ft(:,3));
    ft3 = hypot(ft(:,4), ft(:,5));
    ft_amp = hypot(ft1, ft3);

    TF_dB = 20*log10(max(ft_amp ./ Fw, 1e-300));
end

%% =========================================================
% Helper 3: AFT 立方项批量投影
%% =========================================================
function cubic = cubic_proj_013_batch(U)

    [~, T_mat, T_inv] = get_AFT_matrices_local();

    X_time  = (T_mat * U.').';
    X3_time = X_time.^3;

    cubic = (T_inv * X3_time.').';
end

%% =========================================================
% Helper 4: AFT 矩阵
%% =========================================================
function [N, T_mat, T_inv] = get_AFT_matrices_local()

    persistent pN pT pTinv

    if isempty(pN)
        pN = 64;
        t = (0:pN-1)'*(2*pi/pN);

        dc = ones(pN,1);
        c1 = cos(t);
        s1 = sin(t);
        c3 = cos(3*t);
        s3 = sin(3*t);

        pT = [dc, c1, s1, c3, s3];

        Inv = [dc, 2*c1, 2*s1, 2*c3, 2*s3]';
        pTinv = (1/pN) * Inv;
        pTinv(1,:) = (1/pN) * dc';
    end

    N = pN;
    T_mat = pT;
    T_inv = pTinv;
end