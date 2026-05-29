%% =========================================================
% Validate_04_HBM_vs_ODE_FullModel.m
%
% 目的：
% 验证完整 15 维 HBM-AFT 解是否与原始非线性 ODE 时域积分稳态解一致。
%
% 验证对象：
%   1) xi1(t)
%   2) xi2(t)
%   3) Q(t)
%   4) 基础端传递力 ft(t)
%
% 注意：
%   - 本脚本适用于 lam > 0, kap_e > 0 的完整机电耦合模型。
%   - 如果 HBM 解位于不稳定分支，ODE 长时间积分不会收敛到该 HBM 解。
%   - 因此建议先选稳定、非跳跃区的频率点。
% =========================================================

clear; clc; close all;

%% -----------------------------
% 0) 全局变量
%% -----------------------------
global Fw FixedOmega
Fw = 0.005;
FixedOmega = [];

%% -----------------------------
% 1) 你的常用完整机电参数
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
P.al1 = alpha1 - P.be1;   % 保证 be1 + al1 = alpha1
P.be2 = alpha2;
P.ga1 = gamma1;
P.ga2 = gamma2;
P.mu  = mu;
P.ze1 = 0.05;             % 下层对地阻尼

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

if P.kap_e <= 0
    error('本验证需要 kap_e > 0，因为 ODE 中需要显式求 Q''''。');
end

%% -----------------------------
% 3) 选择验证频率点
%% -----------------------------
% 建议先选几个非跳跃、非强不稳定的频率点。
% 如果某个点 ODE 和 HBM 不一致，不一定是错，可能 HBM 点在不稳定分支。

Omega_list = [0.3, 0.8, 1.2, 2.0];

% ODE 积分参数
nTransientPeriods = 300;   % 积分总周期数
nGridPerPeriod    = 512;   % 最后一个周期采样点数

% Newton 初值
y_guess = zeros(16,1);

results = struct([]);

%% -----------------------------
% 4) 主循环：逐频率验证
%% -----------------------------
for ii = 1:numel(Omega_list)

    Omega = Omega_list(ii);
    fprintf('\n================================================\n');
    fprintf('验证频率 Omega = %.6f\n', Omega);
    fprintf('================================================\n');

    %% 4.1 HBM 求解
    y_guess(16) = Omega;

    [y_hbm, ok, Rn] = newton('nondim_temp2', y_guess, sysP);

    fprintf('HBM Newton ok = %d, residual = %.3e\n', ok, Rn);

    if ~ok || Rn > 1e-6
        warning('Omega=%.6f 处 HBM Newton 残差偏大，跳过该点。', Omega);
        continue;
    end

    % 下一频率用当前解作为初值
    y_guess = y_hbm;

    %% 4.2 HBM 解重构到最后一个周期相位网格
    phi = linspace(0, 2*pi, nGridPerPeriod+1).';
    phi(end) = [];

    hbm = reconstruct_hbm(y_hbm, Omega, phi, sysP);

    %% 4.3 原始 ODE 长时间积分
    T = 2*pi/Omega;
    tEnd   = nTransientPeriods * T;
    tStart = (nTransientPeriods - 1) * T;

    tEval = tStart + phi/Omega;

    % 初始条件：可以用零初值；如果想更快，也可以用 HBM 在 t=0 的状态
    y0_ode = zeros(6,1);
    % 状态顺序：[xi1; xi1dot; xi2; xi2dot; Q; Qdot]

    odeopt = odeset('RelTol',1e-9, ...
                    'AbsTol',1e-10, ...
                    'MaxStep',T/80);

    rhs = @(t, s) fullmodel_ode_rhs(t, s, Omega, Fw, sysP);

    fprintf('开始 ODE 积分：%d periods ...\n', nTransientPeriods);

    sol = ode45(rhs, [0 tEnd], y0_ode, odeopt);

    S = deval(sol, tEval).';

    ode = reconstruct_ode_from_state(S, sysP);

    %% 4.4 误差计算
    e_x1 = rel_rms_error(hbm.x1, ode.x1);
    e_x2 = rel_rms_error(hbm.x2, ode.x2);
    e_Q  = rel_rms_error(hbm.Q,  ode.Q);
    e_ft = rel_rms_error(hbm.ft, ode.ft);

    fprintf('\n--- HBM vs ODE 相对 RMS 误差 ---\n');
    fprintf('e_x1 = %.6e\n', e_x1);
    fprintf('e_x2 = %.6e\n', e_x2);
    fprintf('e_Q  = %.6e\n', e_Q);
    fprintf('e_ft = %.6e\n', e_ft);

    results(ii).Omega = Omega;
    results(ii).e_x1 = e_x1;
    results(ii).e_x2 = e_x2;
    results(ii).e_Q  = e_Q;
    results(ii).e_ft = e_ft;

    %% 4.5 绘图
    figure('Color','w','Position',[100 100 1100 750]);

    subplot(2,2,1);
    plot(phi, hbm.x1, 'k-', 'LineWidth', 2); hold on;
    plot(phi, ode.x1, 'r--', 'LineWidth', 1.5);
    grid on;
    xlabel('\Omega \tau phase');
    ylabel('\xi_1');
    title(sprintf('\\xi_1, \\Omega=%.3f, e=%.2e', Omega, e_x1));
    legend('HBM','ODE','Location','best');

    subplot(2,2,2);
    plot(phi, hbm.x2, 'k-', 'LineWidth', 2); hold on;
    plot(phi, ode.x2, 'b--', 'LineWidth', 1.5);
    grid on;
    xlabel('\Omega \tau phase');
    ylabel('\xi_2');
    title(sprintf('\\xi_2, \\Omega=%.3f, e=%.2e', Omega, e_x2));
    legend('HBM','ODE','Location','best');

    subplot(2,2,3);
    plot(phi, hbm.Q, 'k-', 'LineWidth', 2); hold on;
    plot(phi, ode.Q, 'm--', 'LineWidth', 1.5);
    grid on;
    xlabel('\Omega \tau phase');
    ylabel('Q');
    title(sprintf('Q, \\Omega=%.3f, e=%.2e', Omega, e_Q));
    legend('HBM','ODE','Location','best');

    subplot(2,2,4);
    plot(phi, hbm.ft, 'k-', 'LineWidth', 2); hold on;
    plot(phi, ode.ft, 'g--', 'LineWidth', 1.5);
    grid on;
    xlabel('\Omega \tau phase');
    ylabel('f_t');
    title(sprintf('f_t, \\Omega=%.3f, e=%.2e', Omega, e_ft));
    legend('HBM','ODE','Location','best');

end

%% -----------------------------
% 5) 汇总输出
%% -----------------------------
fprintf('\n\n========== 验证汇总 ==========\n');
fprintf('%10s %14s %14s %14s %14s\n', ...
        'Omega', 'e_x1', 'e_x2', 'e_Q', 'e_ft');

for ii = 1:numel(results)
    if isempty(results(ii).Omega)
        continue;
    end
    fprintf('%10.4f %14.6e %14.6e %14.6e %14.6e\n', ...
            results(ii).Omega, ...
            results(ii).e_x1, ...
            results(ii).e_x2, ...
            results(ii).e_Q, ...
            results(ii).e_ft);
end
fprintf('RMS HBM x1 = %.6e, RMS ODE x1 = %.6e, RMS diff = %.6e\n', ...
        rms(hbm.x1), rms(ode.x1), rms(hbm.x1-ode.x1));
fprintf('RMS HBM ft = %.6e, RMS ODE ft = %.6e, RMS diff = %.6e\n', ...
        rms(hbm.ft), rms(ode.ft), rms(hbm.ft-ode.ft));
fprintf('\n验证完成。\n');

%% =========================================================
% 局部函数 1：HBM 时域重构
%% =========================================================
function out = reconstruct_hbm(y_hbm, Omega, phi, sysP)

    be2 = sysP(2);
    mu  = sysP(3);
    ze1 = sysP(6);
    ga2 = sysP(11);

    x1c = y_hbm(1:5);
    x2c = y_hbm(6:10);
    Qc  = y_hbm(11:15);

    Tmat = basis_matrix(phi);

    D = deriv_matrix(Omega);

    x1dot_c = D * x1c;
    x2dot_c = D * x2c;
    Qdot_c  = D * Qc;

    out.x1 = Tmat * x1c;
    out.x2 = Tmat * x2c;
    out.Q  = Tmat * Qc;

    out.x1dot = Tmat * x1dot_c;
    out.x2dot = Tmat * x2dot_c;
    out.Qdot  = Tmat * Qdot_c;

    out.ft = be2*out.x2 + ga2*out.x2.^3 + 2*mu*ze1*out.x2dot;
end

%% =========================================================
% 局部函数 2：ODE 状态转为输出
%% =========================================================
function out = reconstruct_ode_from_state(S, sysP)

    be2 = sysP(2);
    mu  = sysP(3);
    ze1 = sysP(6);
    ga2 = sysP(11);

    out.x1 = S(:,1);
    out.x1dot = S(:,2);

    out.x2 = S(:,3);
    out.x2dot = S(:,4);

    out.Q = S(:,5);
    out.Qdot = S(:,6);

    out.ft = be2*out.x2 + ga2*out.x2.^3 + 2*mu*ze1*out.x2dot;
end

%% =========================================================
% 局部函数 3：完整模型 ODE 右端
%% =========================================================
function ds = fullmodel_ode_rhs(t, s, Omega, Fw, sysP)

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

    xi1 = s(1);
    v1  = s(2);
    xi2 = s(3);
    v2  = s(4);
    Q   = s(5);
    Qd  = s(6);

    x12 = xi1 - xi2;
    v12 = v1 - v2;

    alpha1 = be1 + al1;

    % 与 nondim_temp2 保持一致：
    % lam > 0 时不启用额外层间阻尼；
    % lam = 0 时可启用断耦合基准阻尼。
    lam_eps = 1e-12;
    zeta12 = 0.0;

    if abs(lam) < lam_eps
        % 如果你的 nondim_temp2 中这里是 0.020412，就保持一致。
        % 如果你改成了自动公式，也在这里同步改。
        zeta12 = 0.020412;
    end

    damp12 = 2*zeta12*v12;

    % 方程 1：
    % xi1'' + alpha1*x12 + ga1*x12^3 + theta*Q' + damp12 = Fw*cos(Omega*t)
    a1 = Fw*cos(Omega*t) ...
         - alpha1*x12 ...
         - ga1*x12^3 ...
         - theta*Qd ...
         - damp12;

    % 方程 2：
    % mu*xi2'' + 2*mu*ze1*xi2' + be2*xi2 + ga2*xi2^3
    % - alpha1*x12 - ga1*x12^3 - theta*Q' - damp12 = 0
    a2 = ( ...
          - 2*mu*ze1*v2 ...
          - be2*xi2 ...
          - ga2*xi2^3 ...
          + alpha1*x12 ...
          + ga1*x12^3 ...
          + theta*Qd ...
          + damp12 ...
         ) / mu;

    % 方程 3：
    % kap_e*Q'' + sigma*Q' + kap_c*Q - theta*(xi1'-xi2') = 0
    Qdd = (theta*v12 - sigma*Qd - kap_c*Q) / kap_e;

    ds = zeros(6,1);
    ds(1) = v1;
    ds(2) = a1;
    ds(3) = v2;
    ds(4) = a2;
    ds(5) = Qd;
    ds(6) = Qdd;
end

%% =========================================================
% 局部函数 4：谐波基矩阵
%% =========================================================
function T = basis_matrix(phi)
    T = [ones(size(phi)), ...
         cos(phi), sin(phi), ...
         cos(3*phi), sin(3*phi)];
end

%% =========================================================
% 局部函数 5：系数域一阶导数矩阵
%% =========================================================
function D = deriv_matrix(Omega)
    D = zeros(5);
    D(2,3) = Omega;
    D(3,2) = -Omega;
    D(4,5) = 3*Omega;
    D(5,4) = -3*Omega;
end

%% =========================================================
% 局部函数 6：相对 RMS 误差
%% =========================================================
function e = rel_rms_error(y1, y2)
    e = rms(y1 - y2) / max(1e-14, rms(y2));
end