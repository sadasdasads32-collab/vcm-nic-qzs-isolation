%% 纯电阻电路 vs 等效电磁阻尼 理论验证脚本
%
% 目的：
%   验证 nondim_temp2_1 中的纯电阻电路支路是否正确。
%
% 验证逻辑：
%   纯电阻电路：
%       sigma*q_dot - theta*(x1_dot - x2_dot) = 0
%
%   可推出：
%       q_dot = theta/sigma * (x1_dot - x2_dot)
%
%   代入机械电磁力：
%       theta*q_dot = theta^2/sigma * (x1_dot - x2_dot)
%
%   因为 theta^2 = lambda，所以纯电阻等效为层间电磁阻尼：
%       c_em = lambda/sigma
%
%   若机械层间阻尼写为：
%       2*zeta_em*(x1_dot - x2_dot)
%
%   则：
%       zeta_em = lambda/(2*sigma)
%
% 本脚本会对多个 Rt 进行对比：
%   1) 电路模型：调用 nondim_temp2_1
%   2) 理论模型：调用本文件内置 nondim_temp2_resistor_equiv
%
% 若两类曲线基本重合，说明纯电阻电路支路正确。
%
% 依赖外部函数：
%   nondim_temp2_1.m
%   newton.m
%   branch_follow2.m

clc; clear; close all;

global Fw FixedOmega
FixedOmega = [];

%% ------------------------------------------------------------
% 1. 物理参数
% ------------------------------------------------------------
Kt = 7.474;
Ke = 7.474;
m1 = 2.2;
k1 = 3000;
R0 = 3.8;

wn = sqrt(k1/m1);

%% ------------------------------------------------------------
% 2. 结构参数
% ------------------------------------------------------------
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

%% ------------------------------------------------------------
% 3. 电路参数：只保留电阻
% ------------------------------------------------------------
% 如果你想复现图中的 lambda = 4.23，可以手动覆盖 P.lam。
% 默认使用物理参数计算：
P.lam = Kt*Ke*wn/(k1*R0);

% 若你想直接验证图中 lambda = 4.23，取消下一行注释：
% P.lam = 4.23;

Rt_list = [3, 6, 9, 12];

P.kap_e = 0.0;
P.kap_c = 0.0;

%% ------------------------------------------------------------
% 4. 求解设置
% ------------------------------------------------------------
Fw = 0.005;

Omega_Start = 10.0;
Omega_Step  = -0.01;
Omega_Next  = Omega_Start + Omega_Step;

n_steps = 3000;

fprintf('\n====================================================\n');
fprintf('纯电阻电路 vs 理论等效电磁阻尼 验证\n');
fprintf('lambda = %.8g\n', P.lam);
fprintf('kap_e = 0, kap_c = 0\n');
fprintf('Rt list = ');
fprintf('%.4g ', Rt_list);
fprintf('\n====================================================\n\n');

%% ------------------------------------------------------------
% 5. 循环计算不同 Rt
% ------------------------------------------------------------
res = struct();

for i = 1:numel(Rt_list)
    Rt = Rt_list(i);
    P.sigma = Rt/R0;

    sysP_circuit = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, ...
                    P.lam, P.kap_e, P.kap_c, P.sigma, P.ga2];

    % 等效阻尼比
    zeta_em = P.lam/(2*P.sigma);

    % 理论等效模型的 sysP：
    % 仍保持 11 个参数格式，但是 kap_e, kap_c, sigma 在等效模型中不使用；
    % lam 和 sigma 用来计算 zeta_em。
    sysP_equiv = sysP_circuit;

    fprintf('Rt = %.6g, sigma = %.6g, c_em = %.6g, zeta_em = %.6g\n', ...
            Rt, P.sigma, P.lam/P.sigma, zeta_em);

    %% -----------------------------
    % A. 纯电阻电路模型：nondim_temp2_1
    % -----------------------------
    y_init = zeros(15,1);
    y_init(end+1) = Omega_Start;

    [x0_full, ok0] = newton('nondim_temp2_1', y_init, sysP_circuit);
    if ~ok0
        error('Rt=%.6g: 电路模型高频起点求解失败。', Rt);
    end
    x0 = x0_full(1:15);

    y_init2 = [x0; Omega_Next];
    [x1_full, ok1] = newton('nondim_temp2_1', y_init2, sysP_circuit);
    if ~ok1
        error('Rt=%.6g: 电路模型第二个点求解失败。', Rt);
    end
    x1 = x1_full(1:15);

    [x_res_circuit, ~] = branch_follow2('nondim_temp2_1', ...
                                        n_steps, ...
                                        Omega_Start, ...
                                        Omega_Next, ...
                                        x0, ...
                                        x1, ...
                                        sysP_circuit);

    [Om_circuit, TF_circuit_dB] = calc_TF_from_xres(x_res_circuit, sysP_circuit, Fw);

    %% -----------------------------
    % B. 理论等效机械阻尼模型
    % -----------------------------
    y_init = zeros(15,1);
    y_init(end+1) = Omega_Start;

    [x0_full, ok0] = newton('nondim_temp2_resistor_equiv', y_init, sysP_equiv);
    if ~ok0
        error('Rt=%.6g: 理论等效模型高频起点求解失败。', Rt);
    end
    x0 = x0_full(1:15);

    y_init2 = [x0; Omega_Next];
    [x1_full, ok1] = newton('nondim_temp2_resistor_equiv', y_init2, sysP_equiv);
    if ~ok1
        error('Rt=%.6g: 理论等效模型第二个点求解失败。', Rt);
    end
    x1 = x1_full(1:15);

    [x_res_equiv, ~] = branch_follow2('nondim_temp2_resistor_equiv', ...
                                      n_steps, ...
                                      Omega_Start, ...
                                      Omega_Next, ...
                                      x0, ...
                                      x1, ...
                                      sysP_equiv);

    [Om_equiv, TF_equiv_dB] = calc_TF_from_xres(x_res_equiv, sysP_equiv, Fw);

    %% -----------------------------
    % C. 插值比较误差
    % -----------------------------
    % 弧长延拓得到的多值曲线不一定是单值函数。
    % 为了给出一个简单数值指标，这里只在唯一频率点上做近似插值比较。
    % 图像重合程度仍然是主要判断依据。
    [Om_c_unique, ia] = unique(round(Om_circuit, 10), 'stable');
    TF_c_unique = TF_circuit_dB(ia);

    [Om_e_unique, ib] = unique(round(Om_equiv, 10), 'stable');
    TF_e_unique = TF_equiv_dB(ib);

    Om_min = max(min(Om_c_unique), min(Om_e_unique));
    Om_max = min(max(Om_c_unique), max(Om_e_unique));

    if Om_min < Om_max
        Om_test = logspace(log10(max(Om_min, 1e-6)), log10(Om_max), 400).';
        TFc_i = interp1(Om_c_unique, TF_c_unique, Om_test, 'linear', NaN);
        TFe_i = interp1(Om_e_unique, TF_e_unique, Om_test, 'linear', NaN);

        ok = isfinite(TFc_i) & isfinite(TFe_i);
        err_dB = TFc_i(ok) - TFe_i(ok);

        max_abs_err = max(abs(err_dB));
        rms_err = sqrt(mean(err_dB.^2));
    else
        max_abs_err = NaN;
        rms_err = NaN;
    end

    fprintf('    dB误差估计: max = %.4g dB, RMS = %.4g dB\n\n', ...
            max_abs_err, rms_err);

    res(i).Rt = Rt;
    res(i).sigma = P.sigma;
    res(i).lambda = P.lam;
    res(i).zeta_em = zeta_em;
    res(i).Om_circuit = Om_circuit;
    res(i).TF_circuit_dB = TF_circuit_dB;
    res(i).Om_equiv = Om_equiv;
    res(i).TF_equiv_dB = TF_equiv_dB;
    res(i).max_abs_err_dB = max_abs_err;
    res(i).rms_err_dB = rms_err;
end

%% ------------------------------------------------------------
% 6. 绘图：每个 Rt 一张对照图
% ------------------------------------------------------------
for i = 1:numel(res)
    figure('Color','w', 'Position',[120 120 760 520]);
    ax = gca;
    hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    set(ax,'XScale','log');

    plot(ax, res(i).Om_circuit, res(i).TF_circuit_dB, 'b-', 'LineWidth', 1.5);
    plot(ax, res(i).Om_equiv,   res(i).TF_equiv_dB,   'r--', 'LineWidth', 1.5);

    yline(ax, 0, 'k--', '0 dB');

    xlabel(ax, '\Omega (log scale)');
    ylabel(ax, 'Force Transmissibility 20log_{10}(|f_t|/f) (dB)');

    title(ax, sprintf('Rt=%.4g: resistor circuit vs equivalent damping, zeta_{em}=%.4g', ...
                      res(i).Rt, res(i).zeta_em));

    legend(ax, ...
           'nondim\_temp2\_1: resistor-only circuit', ...
           'theory: equivalent electromagnetic damping', ...
           'Location', 'best');

    xlim(ax, [0.1, Omega_Start]);
    hold(ax,'off');
end

%% ------------------------------------------------------------
% 7. 绘图：所有 Rt 总览
% ------------------------------------------------------------
figure('Color','w', 'Position',[150 150 900 560]);
ax = gca;
hold(ax,'on'); grid(ax,'on'); box(ax,'on');
set(ax,'XScale','log');

for i = 1:numel(res)
    plot(ax, res(i).Om_circuit, res(i).TF_circuit_dB, '-', 'LineWidth', 1.3, ...
         'DisplayName', sprintf('Circuit Rt=%.4g', res(i).Rt));
    plot(ax, res(i).Om_equiv, res(i).TF_equiv_dB, '--', 'LineWidth', 1.1, ...
         'DisplayName', sprintf('Equiv Rt=%.4g', res(i).Rt));
end

yline(ax, 0, 'k--', '0 dB', 'HandleVisibility','off');

xlabel(ax, '\Omega (log scale)');
ylabel(ax, 'Force Transmissibility 20log_{10}(|f_t|/f) (dB)');
title(ax, sprintf('Validation: resistor-only circuit vs equivalent damping, \\lambda=%.4g', P.lam));
legend(ax, 'Location', 'bestoutside');

xlim(ax, [0.1, Omega_Start]);
hold(ax,'off');

%% ------------------------------------------------------------
% 8. 保存结果
% ------------------------------------------------------------
save('result_resistor_validation.mat', 'res', 'P', 'Rt_list');

fprintf('验证完成，结果已保存到 result_resistor_validation.mat\n');
fprintf('判断标准：同一 Rt 下，蓝色实线和红色虚线应基本重合。\n');

%% ============================================================
% 函数 1：理论等效电磁阻尼模型
% ============================================================
function [z, Jac] = nondim_temp2_resistor_equiv(y, sysP)
    global Fw FixedOmega

    if numel(y) ~= 16
        error('nondim_temp2_resistor_equiv expects y to be 16x1.');
    end
    if numel(sysP) ~= 11
        error('nondim_temp2_resistor_equiv expects sysP to be 11x1.');
    end

    state = y(1:15);

    if isempty(FixedOmega)
        W = y(16);
        current_Fw = Fw;
    else
        W = FixedOmega;
        current_Fw = y(16);
    end

    be1     = sysP(1);
    be2     = sysP(2);
    mu_mass = sysP(3);

    al1 = sysP(4);
    ga1 = sysP(5);
    ze1 = sysP(6);

    lam_phys = sysP(7);
    sigma    = sysP(10);
    ga2      = sysP(11);

    if abs(sigma) < 1e-12
        error('sigma is zero. Equivalent resistor damping requires sigma > 0.');
    end

    % 纯电阻等效层间阻尼：
    % c_em = lambda/sigma = 2*zeta_em
    c_em = lam_phys/sigma;

    x1 = state(1:5);
    x2 = state(6:10);
    q  = state(11:15);

    x12 = x1 - x2;

    cubic12 = cubic_proj_013(x12);
    cubic2  = cubic_proj_013(x2);

    W2 = W^2;

    Mat_Deriv = zeros(5);
    Mat_Deriv(2,3) = W;    Mat_Deriv(3,2) = -W;
    Mat_Deriv(4,5) = 3*W;  Mat_Deriv(5,4) = -3*W;

    Mat_Inertia = diag([0; -W2; -W2; -9*W2; -9*W2]);

    I5 = eye(5);

    x12_dot = Mat_Deriv*x12;
    x2_dot  = Mat_Deriv*x2;

    damp_em = c_em*x12_dot;
    damp2   = (2*mu_mass*ze1)*x2_dot;

    Force_from_upper = (be1+al1)*x12 + ga1*cubic12;

    % 等效机械模型：
    % 原电磁力 theta*q_dot 被替换为 c_em*(x1_dot-x2_dot)
    R1 = Mat_Inertia*x1 ...
         + Force_from_upper ...
         + damp_em;

    R1(2) = R1(2) - current_Fw;

    R2 = mu_mass*(Mat_Inertia*x2) ...
         + damp2 ...
         + be2*x2 + ga2*cubic2 ...
         - Force_from_upper ...
         - damp_em;

    % 为了保持 15 维未知量结构，让 q = 0。
    % 这样不用改 newton 和 branch_follow2。
    R3 = q;

    z = [R1; R2; R3];

    if nargout > 1
        J_cubic_x12 = AFT_GetJac(x12);
        J_cubic_x2  = AFT_GetJac(x2);

        J11 = Mat_Inertia + (be1+al1)*I5 + ga1*J_cubic_x12 + c_em*Mat_Deriv;
        J12 = -(be1+al1)*I5 - ga1*J_cubic_x12 - c_em*Mat_Deriv;
        J13 = zeros(5);

        J21 = -(be1+al1)*I5 - ga1*J_cubic_x12 - c_em*Mat_Deriv;

        J22 = mu_mass*Mat_Inertia ...
              + (2*mu_mass*ze1)*Mat_Deriv ...
              + be2*I5 + ga2*J_cubic_x2 ...
              + (be1+al1)*I5 + ga1*J_cubic_x12 ...
              + c_em*Mat_Deriv;

        J23 = zeros(5);

        J31 = zeros(5);
        J32 = zeros(5);
        J33 = eye(5);

        Jac = [J11, J12, J13;
               J21, J22, J23;
               J31, J32, J33];
    end
end

%% ============================================================
% 函数 2：从 x_res 计算传递率
% ============================================================
function [Om_valid, TF_dB_valid] = calc_TF_from_xres(x_res, sysP, Fw)
    Om  = x_res(16,:).';

    be2 = sysP(2);
    mu  = sysP(3);
    ze1 = sysP(6);
    ga2 = sysP(11);

    x2 = x_res(6:10,:).';

    W = Om;
    x2_dot = zeros(size(x2));
    x2_dot(:,1) = 0;
    x2_dot(:,2) = W .* x2(:,3);
    x2_dot(:,3) = -W .* x2(:,2);
    x2_dot(:,4) = 3*W .* x2(:,5);
    x2_dot(:,5) = -3*W .* x2(:,4);

    x2_cub = cubic_proj_013_batch(x2);

    ft = be2*x2 + ga2*x2_cub + 2*mu*ze1*x2_dot;

    ft1 = hypot(ft(:,2), ft(:,3));
    ft3 = hypot(ft(:,4), ft(:,5));
    ft_amp = hypot(ft1, ft3);

    TF = ft_amp ./ Fw;
    TF_dB = 20*log10(max(TF, 1e-300));

    valid = isfinite(Om) & isfinite(TF_dB) & (Om > 0);
    Om_valid = Om(valid);
    TF_dB_valid = TF_dB(valid);
end

%% ============================================================
% AFT 辅助函数
% ============================================================
function cubic = cubic_proj_013(u)
    [~, T_mat, T_inv] = get_AFT_matrices();
    cubic = T_inv * ((T_mat * u).^3);
end

function J_aft = AFT_GetJac(u)
    [~, T_mat, T_inv] = get_AFT_matrices();
    u_time = T_mat * u;
    df_du = 3 * u_time.^2;
    J_aft = T_inv * (df_du .* T_mat);
end

function cubic = cubic_proj_013_batch(U)
    [~, T_mat, T_inv] = get_AFT_matrices();
    X_time  = (T_mat * U.').';
    X3_time = X_time.^3;
    cubic   = (T_inv * X3_time.').';
end

function [N, T_mat, T_inv] = get_AFT_matrices()
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
