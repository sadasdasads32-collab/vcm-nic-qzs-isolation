%% =========================================================
%  Validate_01_Jacobian_Check.m
%  解析 Jacobian vs 有限差分 Jacobian 的完整 MATLAB 脚本
%
%  目的：
%  验证 nondim_temp2.m 的解析 Jacobian 是否与中心有限差分 Jacobian 一致
%  适用： 
%  判据：
%  rel_err < 1e-6  : 很好
%  rel_err < 1e-5  : 可以接受
%  rel_err > 1e-4  : 大概率 Jacobian 仍有符号或推导错误
% =========================================================

clc; clear; close all;

%% --------------------------------
% 0. 全局变量
%% --------------------------------
global Fw FixedOmega

Fw = 0.005;       % 外激励幅值
FixedOmega = []; % 扫频模式：y(16) = Omega

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

% 组装 sysP
sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, ...
        P.lam, P.kap_e, P.kap_c, P.sigma, P.ga2];

fprintf('\n========== 参数检查 ==========\n');
fprintf('wn      = %.8f rad/s\n', wn);
fprintf('lambda  = %.8f\n', P.lam);
fprintf('kap_e   = %.8f\n', P.kap_e);
fprintf('kap_c   = %.8f\n', P.kap_c);
fprintf('sigma   = %.8f\n', P.sigma);
fprintf('alpha1 = be1 + al1 = %.8f\n', P.be1 + P.al1);
fprintf('alpha2 = be2       = %.8f\n', P.be2);
fprintf('gamma1 = %.8f\n', P.ga1);
fprintf('gamma2 = %.8f\n', P.ga2);

%% --------------------------------
% 3. 构造测试点 y
%% --------------------------------
% 注意：
% 这里不用 Newton 解，而是用一个随机小扰动状态。
% 这样可以避免全零点导致某些非线性 Jacobian 块退化。
% y = [x1(5); x2(5); q(5); Omega]

Omega_test = 0.8;

rng(1);                         % 固定随机种子，便于复现
y = 1e-3 * randn(16,1);          % 小扰动
y(16) = Omega_test;              % 第 16 个变量是 Omega

%% --------------------------------
% 4. 调用 nondim_temp2 得到解析 Jacobian
%% --------------------------------
[R, J_ana] = nondim_temp2(y, sysP);

if any(~isfinite(R))
    error('残差 R 中出现 NaN 或 Inf。');
end

if any(~isfinite(J_ana), 'all')
    error('解析 Jacobian 中出现 NaN 或 Inf。');
end

%% --------------------------------
% 5. 中心有限差分 Jacobian
%% --------------------------------
% nondim_temp2 输出 15 个残差；
% y 有 16 维，但在扫频模式下第 16 维 Omega 作为参数固定，
% Newton 实际只对前 15 个状态变量求导。
% 所以这里有限差分只对 y(1:15) 做。

h = 1e-6;
J_fd = zeros(15,15);

for k = 1:15
    yp = y;
    ym = y;

    yp(k) = yp(k) + h;
    ym(k) = ym(k) - h;

    Rp = nondim_temp2(yp, sysP);
    Rm = nondim_temp2(ym, sysP);

    J_fd(:,k) = (Rp - Rm) / (2*h);
end

%% --------------------------------
% 6. 总体误差
%% --------------------------------
abs_err = norm(J_ana - J_fd, 'fro');
rel_err = abs_err / max(1, norm(J_fd, 'fro'));

max_abs_err = max(abs(J_ana - J_fd), [], 'all');

fprintf('\n========== Jacobian 总体误差 ==========\n');
fprintf('abs_err     = %.6e\n', abs_err);
fprintf('rel_err     = %.6e\n', rel_err);
fprintf('max_abs_err = %.6e\n', max_abs_err);

if rel_err < 1e-6
    fprintf('判定：很好，解析 Jacobian 与有限差分高度一致。\n');
elseif rel_err < 1e-5
    fprintf('判定：可以接受，解析 Jacobian 基本正确。\n');
elseif rel_err < 1e-4
    fprintf('判定：勉强可接受，但建议检查误差最大的块。\n');
else
    fprintf('判定：不通过，Jacobian 很可能仍存在符号或推导错误。\n');
end

%% --------------------------------
% 7. 分块误差诊断
%% --------------------------------
% 状态分块：
% x1 : 1:5
% x2 : 6:10
% q  : 11:15
%
% 残差分块：
% R1 : 1:5
% R2 : 6:10
% R3 : 11:15

blk_names = {'R1-x1','R1-x2','R1-q'; ...
             'R2-x1','R2-x2','R2-q'; ...
             'R3-x1','R3-x2','R3-q'};

row_blocks = {1:5, 6:10, 11:15};
col_blocks = {1:5, 6:10, 11:15};

fprintf('\n========== 分块 Frobenius 相对误差 ==========\n');

for i = 1:3
    for j = 1:3
        rows = row_blocks{i};
        cols = col_blocks{j};

        A = J_ana(rows, cols);
        B = J_fd(rows, cols);

        e_abs = norm(A - B, 'fro');
        e_rel = e_abs / max(1, norm(B, 'fro'));

        fprintf('%6s : abs = %.3e, rel = %.3e\n', ...
                blk_names{i,j}, e_abs, e_rel);
    end
end

%% --------------------------------
% 8. 专门检查 R3 对 x1/x2 的符号
%% --------------------------------
% 这是本次 R3 改符号后最容易出错的地方。
% 如果你已经改成：
%   R3 = ... - theta*D*(x1-x2)
% 那么应有：
%   dR3/dx1 = -theta*D
%   dR3/dx2 = +theta*D

fprintf('\n========== R3-x1 / R3-x2 块误差矩阵 ==========\n');

Err_R3_x1 = J_ana(11:15,1:5) - J_fd(11:15,1:5);
Err_R3_x2 = J_ana(11:15,6:10) - J_fd(11:15,6:10);

disp('Err_R3_x1 = J_ana(11:15,1:5) - J_fd(11:15,1:5)');
disp(Err_R3_x1);

disp('Err_R3_x2 = J_ana(11:15,6:10) - J_fd(11:15,6:10)');
disp(Err_R3_x2);

fprintf('\n========== 验证完成 ==========\n');