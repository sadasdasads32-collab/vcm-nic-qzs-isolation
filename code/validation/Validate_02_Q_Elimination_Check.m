%% =========================================================
%  Validate_02_Q_Elimination_Check.m
%
%  目的：
%  验证完整 15 维模型中求得的 q 谐波系数，
%  是否满足电路方程的频域解析消元关系。
%
%  对应残差符号：
%      R3 = kap_e*q'' + sigma*q' + kap_c*q - theta*(x1'-x2') = 0
%
%  因此频域关系为：
%      (-kap_e*(nOmega)^2 + j*sigma*(nOmega) + kap_c) Q_n
%          = j*(nOmega)*theta*X12_n
%
%      Q_n =
%          j*(nOmega)*theta / (-kap_e*(nOmega)^2 + j*sigma*(nOmega) + kap_c)
%          * X12_n
%
%  判据：
%      基波 Q1 相对误差 < 1e-6  很好
%      三次谐波 Q3 若幅值很小，相对误差可能放大，应同时看绝对误差
% =========================================================

clc; clear; close all;

%% --------------------------------
% 0. 全局变量
%% --------------------------------
global Fw FixedOmega

Fw = 0.005;
FixedOmega = [];   % 扫频模式：y(16) = Omega

%% --------------------------------
% 1. 你的常用物理参数
%% --------------------------------
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

%% --------------------------------
% 2. 机械/QZS 参数
%% --------------------------------
mu   = 0.2;     % 质量比 m2/m1
beta = 2.0;     % 下层竖向线性刚度比
K1   = 1.0;     % 上层水平弹簧刚度比
K2   = 0.5;     % 下层水平弹簧刚度比
U    = 2.0;     % 几何非线性尺度参数
L    = 4/9;     % QZS 长度比

% 由 L = 4/9, K1 = 1 反推 v，使 alpha1 = 0
v = 2.5;

alpha1 = v    - 2*K1*(1-L)/L;
alpha2 = beta - 2*K2*(1-L)/L;
gamma1 = K1/(U^2 * L^3);
gamma2 = K2/(U^2 * L^3);

P.be1 = 1.0;
P.al1 = alpha1 - P.be1;  % 使 be1 + al1 = alpha1
P.be2 = alpha2;
P.ga1 = gamma1;
P.ga2 = gamma2;
P.mu  = mu;
P.ze1 = 0.05;            % 当前代码中对应下层对地阻尼

sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, ...
        P.lam, P.kap_e, P.kap_c, P.sigma, P.ga2];

theta = sqrt(P.lam);

fprintf('\n========== 参数检查 ==========\n');
fprintf('wn      = %.8f rad/s\n', wn);
fprintf('lambda  = %.8f\n', P.lam);
fprintf('theta   = %.8f\n', theta);
fprintf('kap_e   = %.8f\n', P.kap_e);
fprintf('kap_c   = %.8f\n', P.kap_c);
fprintf('sigma   = %.8f\n', P.sigma);
fprintf('alpha1 = %.8f\n', P.be1 + P.al1);
fprintf('alpha2 = %.8f\n', P.be2);
fprintf('gamma1 = %.8f\n', P.ga1);
fprintf('gamma2 = %.8f\n', P.ga2);

%% --------------------------------
% 3. 选择一个验证频率 Omega
%% --------------------------------
% 建议先选一个普通频率点，避开电路共振极近区域。
% 之后可以改成多个 Omega 循环验证。

Omega_test = 0.8;

%% --------------------------------
% 4. 求完整 15 维模型的稳态解
%% --------------------------------
y0 = zeros(16,1);
y0(16) = Omega_test;

[y_sol, ok, Rn] = newton('nondim_temp2', y0, sysP);

fprintf('\n========== Newton 求解 ==========\n');
fprintf('Omega_test = %.8f\n', Omega_test);
fprintf('Newton ok  = %d\n', ok);
fprintf('Residual   = %.6e\n', Rn);

if ~ok
    warning('Newton 未收敛。建议换 Omega_test，或用扫频延拓得到的解来做此验证。');
end

if Rn > 1e-6
    warning('残差偏大，当前点可能不适合用于消元关系验证。');
end
%% --------------------------------
% 5. 提取 x1, x2, q
%% --------------------------------
x1 = y_sol(1:5);
x2 = y_sol(6:10);
q  = y_sol(11:15);

x12 = x1 - x2;
Omega = y_sol(16);

%% --------------------------------
% 6. 定义实系数到复幅值的映射
%% --------------------------------
% 你的谐波展开为：
%   y(t) = y0 + a1 cos(Omega t) + b1 sin(Omega t)
%              + a3 cos(3Omega t) + b3 sin(3Omega t)
%
% 若采用：
%   y(t) = Re{Y1 exp(j Omega t)} + Re{Y3 exp(j 3Omega t)}
%
% 则：
%   Y1 = a1 - j*b1
%   Y3 = a3 - j*b3

X12_1 = x12(2) - 1i*x12(3);
X12_3 = x12(4) - 1i*x12(5);

Q_1 = q(2) - 1i*q(3);
Q_3 = q(4) - 1i*q(5);

%% --------------------------------
% 7. 用电路方程解析消元公式预测 Q1, Q3
%% --------------------------------
n1 = 1;
w1 = n1 * Omega;

D1 = -P.kap_e*w1^2 + 1i*P.sigma*w1 + P.kap_c;
Q_1_pred = (1i*w1*theta / D1) * X12_1;

n3 = 3;
w3 = n3 * Omega;

D3 = -P.kap_e*w3^2 + 1i*P.sigma*w3 + P.kap_c;
Q_3_pred = (1i*w3*theta / D3) * X12_3;

%% --------------------------------
% 8. 误差计算
%% --------------------------------
abs_err_Q1 = abs(Q_1 - Q_1_pred);
rel_err_Q1 = abs_err_Q1 / max(1e-14, abs(Q_1_pred));

abs_err_Q3 = abs(Q_3 - Q_3_pred);
rel_err_Q3 = abs_err_Q3 / max(1e-14, abs(Q_3_pred));

fprintf('\n========== 频域消元关系验证 ==========\n');

fprintf('\n--- 基波 n=1 ---\n');
fprintf('|X12_1|       = %.6e\n', abs(X12_1));
fprintf('|Q_1|         = %.6e\n', abs(Q_1));
fprintf('|Q_1_pred|    = %.6e\n', abs(Q_1_pred));
fprintf('abs_err_Q1   = %.6e\n', abs_err_Q1);
fprintf('rel_err_Q1   = %.6e\n', rel_err_Q1);

fprintf('\n--- 三次谐波 n=3 ---\n');
fprintf('|X12_3|       = %.6e\n', abs(X12_3));
fprintf('|Q_3|         = %.6e\n', abs(Q_3));
fprintf('|Q_3_pred|    = %.6e\n', abs(Q_3_pred));
fprintf('abs_err_Q3   = %.6e\n', abs_err_Q3);
fprintf('rel_err_Q3   = %.6e\n', rel_err_Q3);

%% --------------------------------
% 9. 相位检查
%% --------------------------------
phase_Q1      = angle(Q_1);
phase_Q1_pred = angle(Q_1_pred);
phase_err_Q1  = angle(exp(1i*(phase_Q1 - phase_Q1_pred)));

phase_Q3      = angle(Q_3);
phase_Q3_pred = angle(Q_3_pred);
phase_err_Q3  = angle(exp(1i*(phase_Q3 - phase_Q3_pred)));

fprintf('\n========== 相位误差 ==========\n');
fprintf('phase_err_Q1 = %.6e rad\n', phase_err_Q1);
fprintf('phase_err_Q3 = %.6e rad\n', phase_err_Q3);

%% --------------------------------
% 10. 直接检查 R3 频域残差
%% --------------------------------
% 也就是：
%   Dn*Qn - j*n*Omega*theta*X12n
% 应该接近 0

R3_complex_1 = D1*Q_1 - 1i*w1*theta*X12_1;
R3_complex_3 = D3*Q_3 - 1i*w3*theta*X12_3;

rel_R3_1 = abs(R3_complex_1) / max(1e-14, abs(1i*w1*theta*X12_1));
rel_R3_3 = abs(R3_complex_3) / max(1e-14, abs(1i*w3*theta*X12_3));

fprintf('\n========== 复数形式 R3 残差 ==========\n');
fprintf('abs(R3_complex_1) = %.6e\n', abs(R3_complex_1));
fprintf('rel_R3_1          = %.6e\n', rel_R3_1);
fprintf('abs(R3_complex_3) = %.6e\n', abs(R3_complex_3));
fprintf('rel_R3_3          = %.6e\n', rel_R3_3);

%% --------------------------------
% 11. 判定
%% --------------------------------
fprintf('\n========== 判定 ==========\n');

if rel_err_Q1 < 1e-6 && rel_R3_1 < 1e-6
    fprintf('基波：通过。Q1 与解析消元公式高度一致。\n');
elseif rel_err_Q1 < 1e-4 && rel_R3_1 < 1e-4
    fprintf('基波：基本通过。误差可接受。\n');
else
    fprintf('基波：未通过。请检查 R3 符号、J31/J32、复幅值 a-jb 映射。\n');
end

if abs(Q_3_pred) < 1e-12 && abs(Q_3) < 1e-12
    fprintf('三次谐波：幅值极小，绝对误差更有参考意义；当前可视为通过。\n');
elseif rel_err_Q3 < 1e-6 && rel_R3_3 < 1e-6
    fprintf('三次谐波：通过。Q3 与解析消元公式高度一致。\n');
elseif rel_err_Q3 < 1e-4 && rel_R3_3 < 1e-4
    fprintf('三次谐波：基本通过。误差可接受。\n');
else
    fprintf('三次谐波：未通过。若 |Q3| 极小，请优先看绝对误差；否则检查三次谐波频率 3Omega 的符号。\n');
end

fprintf('\n========== 验证完成 ==========\n');