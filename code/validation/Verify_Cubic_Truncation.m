clear; clc; close all;

%% ============================================================
% HBM vs ODE 硬验证脚本（完整机电一致性版本）
% 目标：
% 1) 用 ODE15s 求稳态周期响应
% 2) 从最后一个周期提取 0/1/3 次谐波系数
% 3) 以该系数为初值求解 HBM
% 4) 比较波形、谐波系数、幅值误差、动力学残差
% =============================================================

%% -----------------------------
% 0. 基础物理参数
% -----------------------------
Kt = 7.474;
Ke = 7.474;
m1 = 2.2;
k1 = 3000;
R0 = 3.8;

wn = sqrt(k1/m1);

Rt  = 2.3674;
Lsh = 0.04065;
Csh = 0.2227;

P = struct();
P.kap_e = Lsh*wn/R0;
P.kap_c = 1/(Csh*R0*wn);
P.sigma = Rt/R0;

%% -----------------------------
% 1. 非线性机械参数
% -----------------------------
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

%% -----------------------------
% 2. 验证测试点
% 建议：
% (A) lam = 0 做退耦对照
% (B) lam > 0 做完整机电验证
% -----------------------------
Test_Omega = 1.2;
Test_Force = 0.0050;

% ===== 这里改验证点 =====
P.lam = 0.05;     % 推荐先用 >0 的值验证完整机电耦合
% P.lam = 0.0;    % 若做退耦对照可改为 0

% 组装 sysP，顺序必须与 nondim_temp2 一致
sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, ...
        P.lam, P.kap_e, P.kap_c, P.sigma, P.ga2];

fprintf('====================================================\n');
fprintf('HBM vs ODE hard validation starts\n');
fprintf('Target: Omega = %.6f, Fw = %.6f\n', Test_Omega, Test_Force);
fprintf('System Params: lam=%.6f, kap_e=%.6f, sigma=%.6f\n', ...
        P.lam, P.kap_e, P.sigma);
fprintf('====================================================\n');

%% -----------------------------
% 3. ODE 时域积分
% -----------------------------
fprintf('[Step 1] Running ODE15s time integration ...\n');

W = Test_Omega;
T = 2*pi/W;

t_end = 1000*T;
y0 = zeros(6,1);

opts = odeset('RelTol',1e-8, ...
              'AbsTol',1e-10, ...
              'MaxStep',T/200);

[t, y_time] = ode15s(@(t,y) sys_ode_current_consistent(t, y, sysP, W, Test_Force), ...
                     [0, t_end], y0, opts);

fprintf('   -> ODE done. Total steps = %d\n', length(t));

%% -----------------------------
% 4. 提取最后一个周期稳态响应
% -----------------------------
fprintf('[Step 2] Extracting last-period steady response ...\n');

idx_1T = t > (t(end) - T);
t_raw  = t(idx_1T);
x1_raw = y_time(idx_1T,1);
x2_raw = y_time(idx_1T,3);
q_raw  = y_time(idx_1T,5);

Ns = 4096;
t_ss = linspace(t_raw(1), t_raw(end), Ns+1).';
t_ss(end) = [];

x1_ss = interp1(t_raw, x1_raw, t_ss, 'pchip');
x2_ss = interp1(t_raw, x2_raw, t_ss, 'pchip');
q_ss  = interp1(t_raw, q_raw , t_ss, 'pchip');

%% -----------------------------
% 5. 提取 ODE 的 0/1/3 次谐波系数
% -----------------------------
fprintf('[Step 3] Computing Fourier coefficients from ODE steady response ...\n');

c_x1 = get_h013_coeff_1T(x1_ss, t_ss, W);
c_x2 = get_h013_coeff_1T(x2_ss, t_ss, W);
c_q  = get_h013_coeff_1T(q_ss , t_ss, W);

A1_x1_ode = hypot(c_x1(2), c_x1(3));
A1_x2_ode = hypot(c_x2(2), c_x2(3));
A1_q_ode  = hypot(c_q(2),  c_q(3));

fprintf('   -> ODE x1 1st-harmonic amplitude = %.8e\n', A1_x1_ode);
fprintf('   -> ODE x2 1st-harmonic amplitude = %.8e\n', A1_x2_ode);
fprintf('   -> ODE q  1st-harmonic amplitude = %.8e\n', A1_q_ode);

%% -----------------------------
% 6. HBM 求解
% -----------------------------
fprintf('[Step 4] Solving HBM algebraic system ...\n');

global FixedOmega Fw
FixedOmega = Test_Omega;
Fw = Test_Force;

y_guess = [c_x1; c_x2; c_q];

fun = @(y15) hbm_residual_15_current(y15, sysP, Test_Force);

options = optimoptions('fsolve', ...
    'Display','iter', ...
    'FunctionTolerance',1e-12, ...
    'StepTolerance',1e-12, ...
    'OptimalityTolerance',1e-12, ...
    'MaxIterations',400, ...
    'MaxFunctionEvaluations',5000);

[y_hbm, fval, exitflag, output] = fsolve(fun, y_guess, options);
Res_Norm = norm(fval);

x1c_hbm = y_hbm(1:5);
x2c_hbm = y_hbm(6:10);
qc_hbm  = y_hbm(11:15);

A1_x1_hbm = hypot(x1c_hbm(2), x1c_hbm(3));
A1_x2_hbm = hypot(x2c_hbm(2), x2c_hbm(3));
A1_q_hbm  = hypot(qc_hbm(2),  qc_hbm(3));

fprintf('   -> HBM exitflag = %d\n', exitflag);
fprintf('   -> HBM iterations = %d\n', output.iterations);
fprintf('   -> Residual norm ||R|| = %.8e\n', Res_Norm);
fprintf('   -> HBM x1 1st-harmonic amplitude = %.8e\n', A1_x1_hbm);
fprintf('   -> HBM x2 1st-harmonic amplitude = %.8e\n', A1_x2_hbm);
fprintf('   -> HBM q  1st-harmonic amplitude = %.8e\n', A1_q_hbm);

%% -----------------------------
% 7. HBM 波形重构
% -----------------------------
fprintf('[Step 5] Reconstructing HBM waveforms ...\n');

phase = W * t_ss;

x1_rec = x1c_hbm(1) ...
       + x1c_hbm(2)*cos(phase) + x1c_hbm(3)*sin(phase) ...
       + x1c_hbm(4)*cos(3*phase) + x1c_hbm(5)*sin(3*phase);

x2_rec = x2c_hbm(1) ...
       + x2c_hbm(2)*cos(phase) + x2c_hbm(3)*sin(phase) ...
       + x2c_hbm(4)*cos(3*phase) + x2c_hbm(5)*sin(3*phase);

q_rec  = qc_hbm(1) ...
       + qc_hbm(2)*cos(phase) + qc_hbm(3)*sin(phase) ...
       + qc_hbm(4)*cos(3*phase) + qc_hbm(5)*sin(3*phase);

%% -----------------------------
% 8. 误差评估
% -----------------------------
fprintf('[Step 6] Evaluating errors ...\n');

Err_x1_A1 = abs(A1_x1_hbm - A1_x1_ode) / max(A1_x1_ode,1e-12) * 100;
Err_x2_A1 = abs(A1_x2_hbm - A1_x2_ode) / max(A1_x2_ode,1e-12) * 100;
Err_q_A1  = abs(A1_q_hbm  - A1_q_ode ) / max(A1_q_ode ,1e-12) * 100;

Err_coeff_x1 = norm(x1c_hbm - c_x1) / max(norm(c_x1),1e-12) * 100;
Err_coeff_x2 = norm(x2c_hbm - c_x2) / max(norm(c_x2),1e-12) * 100;
Err_coeff_q  = norm(qc_hbm  - c_q ) / max(norm(c_q ),1e-12) * 100;

Err_wave_x1 = norm(x1_rec - x1_ss) / max(norm(x1_ss),1e-12) * 100;
Err_wave_x2 = norm(x2_rec - x2_ss) / max(norm(x2_ss),1e-12) * 100;
Err_wave_q  = norm(q_rec  - q_ss ) / max(norm(q_ss ),1e-12) * 100;

fprintf('\n---------------- Error Summary ----------------\n');
fprintf('A1 error x1 = %.6f %%\n', Err_x1_A1);
fprintf('A1 error x2 = %.6f %%\n', Err_x2_A1);
fprintf('A1 error q  = %.6f %%\n', Err_q_A1);
fprintf('Coeff error x1 = %.6f %%\n', Err_coeff_x1);
fprintf('Coeff error x2 = %.6f %%\n', Err_coeff_x2);
fprintf('Coeff error q  = %.6f %%\n', Err_coeff_q);
fprintf('Wave error x1 = %.6f %%\n', Err_wave_x1);
fprintf('Wave error x2 = %.6f %%\n', Err_wave_x2);
fprintf('Wave error q  = %.6f %%\n', Err_wave_q);

%% -----------------------------
% 9. 用 HBM 重构波形代回 ODE，检查动力学残差
% -----------------------------
fprintf('[Step 7] Checking ODE residual from HBM reconstruction ...\n');

[x1h, x1dh, x1ddh, x2h, x2dh, x2ddh, qh, qdh, qddh] = ...
    reconstruct_h013_with_derivatives(x1c_hbm, x2c_hbm, qc_hbm, W, t_ss);

[r1, r2, r3] = ode_residual_from_reconstruction( ...
    x1h, x1dh, x1ddh, x2h, x2dh, x2ddh, qh, qdh, qddh, ...
    sysP, W, Test_Force, t_ss);

RMS_r1 = sqrt(mean(r1.^2));
RMS_r2 = sqrt(mean(r2.^2));
RMS_r3 = sqrt(mean(r3.^2));

fprintf('ODE residual RMS r1 = %.8e\n', RMS_r1);
fprintf('ODE residual RMS r2 = %.8e\n', RMS_r2);
fprintf('ODE residual RMS r3 = %.8e\n', RMS_r3);

%% -----------------------------
% 10. 判定
% -----------------------------
fprintf('\n================ Final Verdict ================\n');

if Res_Norm < 1e-8 && max([Err_wave_x1, Err_wave_x2, Err_wave_q]) < 3
    fprintf('STATUS: [GOOD MATCH] HBM-AFT framework is numerically validated.\n');
elseif Res_Norm < 1e-6 && max([Err_wave_x1, Err_wave_x2, Err_wave_q]) < 8
    fprintf('STATUS: [ACCEPTABLE] Main trend is captured, but higher harmonics may be insufficient.\n');
else
    fprintf('STATUS: [CHECK NEEDED] Try larger harmonic order or re-check nondim_temp2 consistency.\n');
end

fprintf('====================================================\n');

%% -----------------------------
% 11. 作图
% -----------------------------
figure(200); clf; set(gcf,'Color','w','Position',[120,80,1000,780]);

subplot(3,1,1);
plot(t_ss, x1_ss, 'b-', 'LineWidth',3); hold on;
plot(t_ss, x1_rec, 'r--', 'LineWidth',1.8);
xlim([t_ss(1), t_ss(1)+3*T]);
xlabel('Time \tau'); ylabel('\xi_1');
title(sprintf('x_1: ODE vs HBM (\\Omega=%.3f, F_w=%.4f)', Test_Omega, Test_Force));
legend('ODE15s','HBM (0,1,3)','Location','best');
grid on; set(gca,'FontSize',12);

subplot(3,1,2);
plot(t_ss, x2_ss, 'b-', 'LineWidth',3); hold on;
plot(t_ss, x2_rec, 'r--', 'LineWidth',1.8);
xlim([t_ss(1), t_ss(1)+3*T]);
xlabel('Time \tau'); ylabel('\xi_2');
title('x_2: ODE vs HBM reconstruction');
legend('ODE15s','HBM (0,1,3)','Location','best');
grid on; set(gca,'FontSize',12);

subplot(3,1,3);
plot(t_ss, q_ss, 'b-', 'LineWidth',3); hold on;
plot(t_ss, q_rec, 'r--', 'LineWidth',1.8);
xlim([t_ss(1), t_ss(1)+3*T]);
xlabel('Time \tau'); ylabel('q');
title('q: ODE vs HBM reconstruction');
legend('ODE15s','HBM (0,1,3)','Location','best');
grid on; set(gca,'FontSize',12);

figure(201); clf; set(gcf,'Color','w','Position',[180,120,900,420]);
bar([A1_x1_hbm, A1_x1_ode; A1_x2_hbm, A1_x2_ode; A1_q_hbm, A1_q_ode]);
set(gca,'XTickLabel',{'x_1','x_2','q'},'FontSize',12);
ylabel('1st harmonic amplitude');
legend('HBM','ODE','Location','best');
title(sprintf(['First harmonic comparison: err_{x1}=%.3f%%, err_{x2}=%.3f%%, err_q=%.3f%%'], ...
      Err_x1_A1, Err_x2_A1, Err_q_A1));
grid on;

figure(202); clf; set(gcf,'Color','w','Position',[200,160,900,420]);
bar([Err_coeff_x1, Err_coeff_x2, Err_coeff_q; ...
     Err_wave_x1,  Err_wave_x2,  Err_wave_q]');
set(gca,'XTickLabel',{'x_1','x_2','q'},'FontSize',12);
ylabel('Relative error (%)');
legend('Coefficient error','Waveform error','Location','best');
title('HBM vs ODE errors');
grid on;

figure(203); clf; set(gcf,'Color','w','Position',[220,180,1000,700]);

subplot(3,1,1);
plot(t_ss, r1, 'k-', 'LineWidth',1.5);
xlim([t_ss(1), t_ss(1)+3*T]);
ylabel('r_1');
title(sprintf('ODE residual from HBM reconstruction, RMS = %.3e', RMS_r1));
grid on;

subplot(3,1,2);
plot(t_ss, r2, 'k-', 'LineWidth',1.5);
xlim([t_ss(1), t_ss(1)+3*T]);
ylabel('r_2');
title(sprintf('ODE residual from HBM reconstruction, RMS = %.3e', RMS_r2));
grid on;

subplot(3,1,3);
plot(t_ss, r3, 'k-', 'LineWidth',1.5);
xlim([t_ss(1), t_ss(1)+3*T]);
xlabel('Time \tau'); ylabel('r_3');
title(sprintf('ODE residual from HBM reconstruction, RMS = %.3e', RMS_r3));
grid on;

%% ============================================================
% 局部函数 1：HBM 残差包装
% ============================================================
function R15 = hbm_residual_15_current(y15, sysP, Fw_fixed)
    y16 = [y15(:); Fw_fixed];
    R15 = nondim_temp2(y16, sysP);
end

%% ============================================================
% 局部函数 2：与 nondim_temp2 一致的 ODE 方程
% 注意：这里显式包含 lam≈0 时 zeta12=0.05 的逻辑
% ============================================================
function dydt = sys_ode_current_consistent(t, y, sysP, W, Fw)

    be1   = sysP(1);
    be2   = sysP(2);
    mu    = sysP(3);
    al1   = sysP(4);
    ga1   = sysP(5);
    ze1   = sysP(6);
    lam   = sysP(7);
    kap_e = sysP(8);
    kap_c = sysP(9);
    sigma = sysP(10);
    ga2   = sysP(11);

    x1  = y(1); x1d = y(2);
    x2  = y(3); x2d = y(4);
    q   = y(5); qd  = y(6);

    theta = sqrt(max(0, lam));

    lam_eps = 1e-12;
    zeta12 = 0.0;
    if abs(lam) < lam_eps
        zeta12 = 0.05;
    end

    x12  = x1 - x2;
    x12d = x1d - x2d;

    dydt = zeros(6,1);

    % 上层
    dydt(1) = x1d;
    dydt(2) = Fw*cos(W*t) ...
            - (be1 + al1)*x12 ...
            - ga1*x12^3 ...
            - theta*qd ...
            - 2*zeta12*x12d;

    % 下层
    dydt(3) = x2d;
    dydt(4) = ((be1 + al1)*x12 ...
             + ga1*x12^3 ...
             + theta*qd ...
             + 2*zeta12*x12d ...
             - 2*mu*ze1*x2d ...
             - be2*x2 ...
             - ga2*x2^3) / mu;

    % 电路
    dydt(5) = qd;
    dydt(6) = (-sigma*qd - kap_c*q - theta*x12d) / kap_e;
end

%% ============================================================
% 局部函数 3：最后一个周期上提取 0/1/3 次谐波
% ============================================================
function coeff = get_h013_coeff_1T(u, t, W)
    T = 2*pi/W;
    a0 = 1/T * trapz(t, u);
    a1 = 2/T * trapz(t, u .* cos(W*t));
    b1 = 2/T * trapz(t, u .* sin(W*t));
    a3 = 2/T * trapz(t, u .* cos(3*W*t));
    b3 = 2/T * trapz(t, u .* sin(3*W*t));
    coeff = [a0; a1; b1; a3; b3];
end

%% ============================================================
% 局部函数 4：重构位移/电荷及其一二阶导数
% coeff = [a0; a1; b1; a3; b3]
% ============================================================
function [x1, x1d, x1dd, x2, x2d, x2dd, q, qd, qdd] = ...
    reconstruct_h013_with_derivatives(c1, c2, cq, W, t)

    wt = W*t;

    x1 = c1(1) + c1(2)*cos(wt) + c1(3)*sin(wt) + c1(4)*cos(3*wt) + c1(5)*sin(3*wt);
    x2 = c2(1) + c2(2)*cos(wt) + c2(3)*sin(wt) + c2(4)*cos(3*wt) + c2(5)*sin(3*wt);
    q  = cq(1) + cq(2)*cos(wt) + cq(3)*sin(wt) + cq(4)*cos(3*wt) + cq(5)*sin(3*wt);

    x1d = -W*c1(2)*sin(wt) + W*c1(3)*cos(wt) - 3*W*c1(4)*sin(3*wt) + 3*W*c1(5)*cos(3*wt);
    x2d = -W*c2(2)*sin(wt) + W*c2(3)*cos(wt) - 3*W*c2(4)*sin(3*wt) + 3*W*c2(5)*cos(3*wt);
    qd  = -W*cq(2)*sin(wt) + W*cq(3)*cos(wt) - 3*W*cq(4)*sin(3*wt) + 3*W*cq(5)*cos(3*wt);

    x1dd = -W^2*c1(2)*cos(wt) - W^2*c1(3)*sin(wt) - 9*W^2*c1(4)*cos(3*wt) - 9*W^2*c1(5)*sin(3*wt);
    x2dd = -W^2*c2(2)*cos(wt) - W^2*c2(3)*sin(wt) - 9*W^2*c2(4)*cos(3*wt) - 9*W^2*c2(5)*sin(3*wt);
    qdd  = -W^2*cq(2)*cos(wt) - W^2*cq(3)*sin(wt) - 9*W^2*cq(4)*cos(3*wt) - 9*W^2*cq(5)*sin(3*wt);
end

%% ============================================================
% 局部函数 5：HBM 重构解代回 ODE 的残差
% ============================================================
function [r1, r2, r3] = ode_residual_from_reconstruction( ...
    x1, x1d, x1dd, x2, x2d, x2dd, q, qd, qdd, sysP, W, Fw, t)

    be1   = sysP(1);
    be2   = sysP(2);
    mu    = sysP(3);
    al1   = sysP(4);
    ga1   = sysP(5);
    ze1   = sysP(6);
    lam   = sysP(7);
    kap_e = sysP(8);
    kap_c = sysP(9);
    sigma = sysP(10);
    ga2   = sysP(11);

    theta = sqrt(max(0, lam));

    lam_eps = 1e-12;
    zeta12 = 0.0;
    if abs(lam) < lam_eps
        zeta12 = 0.05;
    end

    x12  = x1 - x2;
    x12d = x1d - x2d;

    r1 = x1dd + (be1 + al1)*x12 + ga1*x12.^3 + theta*qd + 2*zeta12*x12d - Fw*cos(W*t);
    r2 = mu*x2dd + 2*mu*ze1*x2d + be2*x2 + ga2*x2.^3 - (be1 + al1)*x12 - ga1*x12.^3 - theta*qd - 2*zeta12*x12d;
    r3 = kap_e*qdd + sigma*qd + kap_c*q + theta*x12d;
end