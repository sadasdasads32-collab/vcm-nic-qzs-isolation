%% Validate_HBM_vs_ODE45.m
% =========================================================
% 目的：
% 验证 HBM-AFT 求得的完整 15 维周期解，是否与原始 ODE45
% 时域积分结果一致。
%
% 修正点：
% ODE 电路方程与当前 nondim_temp2 保持一致：
%   kap_e*Q'' + sigma*Q' + kap_c*Q - theta*(x1'-x2') = 0
%
% 即：
%   Q'' = (theta*(x1'-x2') - sigma*Q' - kap_c*Q)/kap_e
%
% 状态顺序：
%   S = [x1; x2; Q; v1; v2; Qdot]
% =========================================================

clc; clear; close all;
init_path();

%% -----------------------------
% 1) 常用完整机电参数
%% -----------------------------
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

%% -----------------------------
% 2) 机械/QZS 参数
%% -----------------------------
mu   = 0.2;
beta = 2.0;
K1   = 1.0;
K2   = 0.5;
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
global Fw FixedOmega
Fw = 0.005;
FixedOmega = [];

Omega_test = 1.2;

fprintf('\n==== HBM vs ODE45 验证开始 ====\n');
fprintf('Omega = %.6f, Fw = %.6f\n', Omega_test, Fw);
fprintf('alpha1 = be1 + al1 = %.6f\n', P.be1 + P.al1);
fprintf('lambda = %.6f, theta = %.6f\n', P.lam, sqrt(P.lam));
fprintf('kap_e = %.6f, kap_c = %.6f, sigma = %.6f\n', ...
        P.kap_e, P.kap_c, P.sigma);

if P.kap_e <= 0
    error('本脚本需要 kap_e > 0，因为 ODE 中需要显式计算 Q''''。');
end

%% ---------------------------------------------------------
% 2) HBM 求解
%% ---------------------------------------------------------
y_init = zeros(16,1);
y_init(16) = Omega_test;

[y_hbm, ok, Rn] = newton('nondim_temp2', y_init, sysP);

fprintf('\n---- HBM Newton ----\n');
fprintf('ok = %d, residual = %.6e\n', ok, Rn);

if ~ok || Rn > 1e-6
    warning('HBM Newton 残差偏大，后续对比可能不可靠。');
end

x1_hbm = y_hbm(1:5);
x2_hbm = y_hbm(6:10);
Q_hbm  = y_hbm(11:15);

%% ---------------------------------------------------------
% 3) 由 HBM 解构造 ODE 初始条件
%% ---------------------------------------------------------
% 谐波形式：
% x(t) = dc + a1*cos(Wt) + b1*sin(Wt)
%      + a3*cos(3Wt) + b3*sin(3Wt)
%
% t=0:
% x(0)    = dc + a1 + a3
% xdot(0) = W*b1 + 3W*b3

get_IC = @(C, W) [C(1) + C(2) + C(4); ...
                  W*C(3) + 3*W*C(5)];

IC_x1 = get_IC(x1_hbm, Omega_test);
IC_x2 = get_IC(x2_hbm, Omega_test);
IC_Q  = get_IC(Q_hbm,  Omega_test);

% 状态顺序：
% S = [x1; x2; Q; v1; v2; Qdot]
S0 = [IC_x1(1);
      IC_x2(1);
      IC_Q(1);
      IC_x1(2);
      IC_x2(2);
      IC_Q(2)];

%% ---------------------------------------------------------
% 4) ODE45 积分
%% ---------------------------------------------------------
T_period = 2*pi/Omega_test;

num_periods = 100;
tEnd = num_periods*T_period;

% 为了对比更干净，固定取最后 1 个周期的均匀采样点
nGrid = 1000;
t_last_start = (num_periods - 1)*T_period;
t_eval = linspace(t_last_start, tEnd, nGrid).';

odeopt = odeset('RelTol',1e-9, ...
                'AbsTol',1e-10, ...
                'MaxStep',T_period/100);

fprintf('\n---- ODE45 积分 ----\n');
fprintf('积分周期数 = %d\n', num_periods);

rhs = @(t, S) sys_ode_fixed(t, S, sysP, Omega_test, Fw);

sol = ode45(rhs, [0, tEnd], S0, odeopt);
S_ode = deval(sol, t_eval).';

%% ---------------------------------------------------------
% 5) 在相同时间点重构 HBM 时域响应
%% ---------------------------------------------------------
t = t_eval;
W = Omega_test;

x1_hbm_time = reconstruct_013(x1_hbm, W, t);
x2_hbm_time = reconstruct_013(x2_hbm, W, t);
Q_hbm_time  = reconstruct_013(Q_hbm,  W, t);

x1_ode = S_ode(:,1);
x2_ode = S_ode(:,2);
Q_ode  = S_ode(:,3);

%% ---------------------------------------------------------
% 6) 计算误差
%% ---------------------------------------------------------
err_x1_max = max(abs(x1_ode - x1_hbm_time));
err_x2_max = max(abs(x2_ode - x2_hbm_time));
err_Q_max  = max(abs(Q_ode  - Q_hbm_time));

err_x1_rms = rms(x1_ode - x1_hbm_time) / max(1e-14, rms(x1_ode));
err_x2_rms = rms(x2_ode - x2_hbm_time) / max(1e-14, rms(x2_ode));
err_Q_rms  = rms(Q_ode  - Q_hbm_time)  / max(1e-14, rms(Q_ode));

fprintf('\n==== 验证结果 ====\n');
fprintf('x1 最大绝对误差      = %.6e\n', err_x1_max);
fprintf('x2 最大绝对误差      = %.6e\n', err_x2_max);
fprintf('Q  最大绝对误差      = %.6e\n', err_Q_max);

fprintf('x1 相对 RMS 误差     = %.6e\n', err_x1_rms);
fprintf('x2 相对 RMS 误差     = %.6e\n', err_x2_rms);
fprintf('Q  相对 RMS 误差     = %.6e\n', err_Q_rms);

if max([err_x1_rms, err_x2_rms, err_Q_rms]) < 1e-3
    fprintf('--> 通过：HBM 与 ODE45 稳态响应高度一致。\n');
elseif max([err_x1_rms, err_x2_rms, err_Q_rms]) < 1e-2
    fprintf('--> 基本通过：误差在可接受范围内。\n');
else
    fprintf('--> 误差偏大：请检查该频点是否位于不稳定分支、跳跃区，或谐波截断是否不足。\n');
end

%% ---------------------------------------------------------
% 7) 绘图
%% ---------------------------------------------------------
phase = mod(W*t, 2*pi);

figure('Name','HBM vs ODE45 Validation FIXED','Color','w', ...
       'Position',[100 100 1050 700]);

subplot(3,1,1);
plot(phase, x1_ode, 'k-', 'LineWidth', 2); hold on;
plot(phase, x1_hbm_time, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('\Omega t');
ylabel('x_1');
title(sprintf('x_1 comparison, rel RMS = %.2e', err_x1_rms));
legend('ODE45','HBM','Location','best');

subplot(3,1,2);
plot(phase, x2_ode, 'k-', 'LineWidth', 2); hold on;
plot(phase, x2_hbm_time, 'b--', 'LineWidth', 1.5);
grid on;
xlabel('\Omega t');
ylabel('x_2');
title(sprintf('x_2 comparison, rel RMS = %.2e', err_x2_rms));
legend('ODE45','HBM','Location','best');

subplot(3,1,3);
plot(phase, Q_ode, 'k-', 'LineWidth', 2); hold on;
plot(phase, Q_hbm_time, 'g--', 'LineWidth', 1.5);
grid on;
xlabel('\Omega t');
ylabel('Q');
title(sprintf('Q comparison, rel RMS = %.2e', err_Q_rms));
legend('ODE45','HBM','Location','best');

fprintf('\n==== 验证完成 ====\n');

%% =========================================================
% 局部函数 1：重构 0/1/3 谐波时域响应
%% =========================================================
function x = reconstruct_013(C, W, t)
    x = C(1) ...
      + C(2)*cos(W*t) + C(3)*sin(W*t) ...
      + C(4)*cos(3*W*t) + C(5)*sin(3*W*t);
end

%% =========================================================
% 局部函数 2：修正后的原始 ODE
%% =========================================================
function dS = sys_ode_fixed(t, S, sysP, W, Fw)

    be1 = sysP(1);
    be2 = sysP(2);
    mu  = sysP(3);

    al1 = sysP(4);
    ga1 = sysP(5);
    ze1 = sysP(6);

    lam   = sysP(7);
    kap_e = sysP(8);
    kap_c = sysP(9);
    sigma = sysP(10);
    ga2   = sysP(11);

    theta = sqrt(max(0, lam));

    % 状态顺序：
    % S = [x1; x2; Q; v1; v2; Qdot]
    x1 = S(1);
    x2 = S(2);
    Q  = S(3);

    v1 = S(4);
    v2 = S(5);
    Qd = S(6);

    x12 = x1 - x2;
    v12 = v1 - v2;

    alpha1 = be1 + al1;

    Force_upper = alpha1*x12 + ga1*x12^3;
    Force_lower = be2*x2 + ga2*x2^3;

    force_em = theta * Qd;

    % 与 nondim_temp2 保持一致：
    % lam > 0 时无额外层间阻尼；
    % lam = 0 时才启用断耦合基准阻尼。
    zeta12 = 0.0;
    if abs(lam) < 1e-12
        zeta12 = 0.05;
    end

    damp12 = 2*zeta12*v12;
    damp2  = 2*mu*ze1*v2;

    % 机械方程 1：
    % x1'' + Force_upper + theta*Q' + damp12 = Fw*cos(Wt)
    a1 = Fw*cos(W*t) ...
       - Force_upper ...
       - force_em ...
       - damp12;

    % 机械方程 2：
    % mu*x2'' + damp2 + Force_lower - Force_upper - theta*Q' - damp12 = 0
    a2 = (-damp2 ...
          - Force_lower ...
          + Force_upper ...
          + force_em ...
          + damp12) / mu;

    % 电路方程，修正后与当前 nondim_temp2 一致：
    % kap_e*Q'' + sigma*Q' + kap_c*Q - theta*(x1'-x2') = 0
    Qdd = (theta*v12 - sigma*Qd - kap_c*Q) / kap_e;

    dS = zeros(6,1);
    dS(1) = v1;
    dS(2) = v2;
    dS(3) = Qd;
    dS(4) = a1;
    dS(5) = a2;
    dS(6) = Qdd;
end