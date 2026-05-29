%% =========================================================
% Plot_Chapter4_1_ComplexOperator.m
%
% 目的：
% 绘制第 4.1 节复动力学算子 K(Omega) 的实部、虚部及物理分区
%
% 当前符号约定：
% 电路方程：
%   kap_e*Q'' + sigma*Q' + kap_c*Q - theta*(x1' - x2') = 0
%
% 频域消元后：
%   F_em = K(Omega) * (xi1 - xi2)
%
%   K(Omega) =
%       theta^2 * Omega^2
%       --------------------------------
%       kap_e*Omega^2 - j*sigma*Omega - kap_c
%
% 写成：
%   K = K_r + j K_i
%
% 其中：
%   K_r = theta^2*Omega^2*(kap_e*Omega^2 - kap_c)
%         / [ (kap_e*Omega^2 - kap_c)^2 + (sigma*Omega)^2 ]
%
%   K_i = theta^2*sigma*Omega^3
%         / [ (kap_e*Omega^2 - kap_c)^2 + (sigma*Omega)^2 ]
% =========================================================

clear; clc; close all;

%% ---------------------------------------------------------
% 1) 使用你前面验证中的常用物理参数
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

lam   = Kt*Ke*wn/(k1*R0);
kap_e = Lsh*wn/R0;
kap_c = 1/(Csh*R0*wn);
sigma = Rt/R0;

theta = sqrt(lam);

fprintf('\n========== Operator parameters ==========\n');
fprintf('lambda  = %.8f\n', lam);
fprintf('theta   = %.8f\n', theta);
fprintf('kap_e   = %.8f\n', kap_e);
fprintf('kap_c   = %.8f\n', kap_c);
fprintf('sigma   = %.8f\n', sigma);

%% ---------------------------------------------------------
% 2) 频率范围
%% ---------------------------------------------------------
Omega_min = 0.05;
Omega_max = 10.0;
N = 3000;

Omega = logspace(log10(Omega_min), log10(Omega_max), N).';

%% ---------------------------------------------------------
% 3) 复算子 K(Omega)
%% ---------------------------------------------------------
Delta = kap_e*Omega.^2 - kap_c;
Den   = Delta.^2 + (sigma*Omega).^2;

Kr = theta^2 .* Omega.^2 .* Delta ./ Den;
Ki = theta^2 .* sigma .* Omega.^3 ./ Den;

K_complex = Kr + 1i*Ki;

% 电路特征频率
Omega_e = sqrt(kap_c/kap_e);

% 等效阻尼和等效惯性/刚度解释量
c_eq = Ki ./ Omega;               % 等效阻尼
m_eq = -Kr ./ (Omega.^2);          % 若解释为惯性项 -Omega^2*m_eq
k_eq = Kr;                         % 若解释为位移同相刚度项

fprintf('Omega_e = sqrt(kap_c/kap_e) = %.8f\n', Omega_e);

%% ---------------------------------------------------------
% 4) 图 1：K_r 与 K_i
%% ---------------------------------------------------------
figure('Color','w','Position',[120 120 900 520]);
hold on; grid on; box on;
set(gca,'XScale','log');

plot(Omega, Kr, 'b-', 'LineWidth', 2.0);
plot(Omega, Ki, 'r--', 'LineWidth', 2.0);

xline(Omega_e, 'k-.', 'LineWidth', 1.5);

yline(0, 'k-', 'LineWidth', 0.8);

xlabel('\Omega');
ylabel('K(\Omega)');
title('Real and imaginary parts of the complex dynamic operator');

legend('K_r(\Omega): real part', ...
       'K_i(\Omega): imaginary part', ...
       '\Omega_e = \surd(\kappa_c/\kappa_e)', ...
       'Location','best');

% 添加文字标注
yl = ylim;
text(Omega_e*1.08, yl(2)*0.82, ...
    sprintf('\\Omega_e = %.3f', Omega_e), ...
    'FontSize', 11, 'Color','k');

text(Omega_min*1.25, yl(1)*0.65 + yl(2)*0.35, ...
    'virtual inertia region', ...
    'FontSize', 10, 'Color','b');

text(Omega_e*1.15, yl(1)*0.25 + yl(2)*0.75, ...
    'damping-shaping region', ...
    'FontSize', 10, 'Color','r');

text(Omega_e*3.0, yl(1)*0.65 + yl(2)*0.35, ...
    'stiffness-like region', ...
    'FontSize', 10, 'Color','b');

xlim([Omega_min Omega_max]);

%% ---------------------------------------------------------
% 5) 图 2：三种等效解释量 m_eq, c_eq, k_eq
%% ---------------------------------------------------------
figure('Color','w','Position',[120 120 1000 760]);

subplot(3,1,1);
hold on; grid on; box on;
set(gca,'XScale','log');
plot(Omega, m_eq, 'b-', 'LineWidth', 1.8);
xline(Omega_e, 'k-.', 'LineWidth', 1.2);
xlabel('\Omega');
ylabel('m_{eq}(\Omega)');
title('Equivalent virtual inertia interpretation: m_{eq} = -K_r/\Omega^2');

subplot(3,1,2);
hold on; grid on; box on;
set(gca,'XScale','log');
plot(Omega, c_eq, 'r-', 'LineWidth', 1.8);
xline(Omega_e, 'k-.', 'LineWidth', 1.2);
xlabel('\Omega');
ylabel('c_{eq}(\Omega)');
title('Equivalent damping interpretation: c_{eq} = K_i/\Omega');

subplot(3,1,3);
hold on; grid on; box on;
set(gca,'XScale','log');
plot(Omega, k_eq, 'm-', 'LineWidth', 1.8);
xline(Omega_e, 'k-.', 'LineWidth', 1.2);
yline(0, 'k-', 'LineWidth', 0.8);
xlabel('\Omega');
ylabel('k_{eq}(\Omega)');
title('Equivalent stiffness interpretation: k_{eq} = K_r');

xlim([Omega_min Omega_max]);

%% ---------------------------------------------------------
% 6) 图 3：归一化后的 K_r, K_i，便于论文展示形状
%% ---------------------------------------------------------
Kr_norm = Kr ./ max(abs(Kr));
Ki_norm = Ki ./ max(abs(Ki));

figure('Color','w','Position',[120 120 900 500]);
hold on; grid on; box on;
set(gca,'XScale','log');

plot(Omega, Kr_norm, 'b-', 'LineWidth', 2.0);
plot(Omega, Ki_norm, 'r--', 'LineWidth', 2.0);
xline(Omega_e, 'k-.', 'LineWidth', 1.5);
yline(0, 'k-', 'LineWidth', 0.8);

xlabel('\Omega');
ylabel('Normalized value');
title('Normalized real and imaginary parts of K(\Omega)');
legend('Normalized K_r', 'Normalized K_i', ...
       '\Omega_e', 'Location','best');

xlim([Omega_min Omega_max]);

%% ---------------------------------------------------------
% 7) 输出关键极限值，方便写论文
%% ---------------------------------------------------------
Kr_low_approx_coeff = -theta^2/kap_c;
m_eq_low = theta^2/kap_c;
Kr_high_limit = theta^2/kap_e;

fprintf('\n========== Asymptotic interpretation ==========\n');
fprintf('Low-frequency: K_r ≈ -(theta^2/kap_c)*Omega^2\n');
fprintf('theta^2/kap_c = %.8f\n', theta^2/kap_c);
fprintf('So virtual inertia m_eq_low ≈ %.8f\n', m_eq_low);
fprintf('High-frequency: K_r -> theta^2/kap_e = %.8f\n', Kr_high_limit);

%% ---------------------------------------------------------
% 8) 导出数据到工作区
%% ---------------------------------------------------------
OperatorData = table(Omega, Kr, Ki, c_eq, m_eq, k_eq);
assignin('base', 'OperatorData', OperatorData);
assignin('base', 'Omega_e', Omega_e);
assignin('base', 'Kr', Kr);
assignin('base', 'Ki', Ki);
assignin('base', 'c_eq', c_eq);
assignin('base', 'm_eq', m_eq);
assignin('base', 'k_eq', k_eq);

fprintf('\nDone. Operator data exported to workspace as OperatorData.\n');