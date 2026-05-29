%% =========================================================
% Validate_05_BranchFollow_Check.m
%
% 目的：
% 验证 branch_follow2 得到的 FRF 分支是否由真实 HBM 解组成，
% 并检查分支点用 fixed-Omega Newton 重新校正后是否保持一致。
%
% 验证内容：
% 1) 分支覆盖范围
% 2) 所有延拓点残差 norm(R)
% 3) 抽样点 branch 解 vs 定频 Newton 校正解
% 4) 抽样点传递力 TF 差异
%
% 依赖：
%   nondim_temp2.m
%   newton.m
%   branch_follow2.m
%   branch_aux2.m
% =========================================================

clear; clc; close all;

%% -----------------------------
% 0) 全局变量
%% -----------------------------
global Fw FixedOmega ParamMin ParamMax

Fw = 0.005;
FixedOmega = [];     % 扫频模式：y(16)=Omega

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

%% -----------------------------
% 3) 延拓设置
%% -----------------------------
% 这里建议先从高频向低频扫，通常高频点更容易收敛。
Omega_Start = 10.0;
Omega_End   = 0.2;
Omega_Step  = -0.01;

nStepsArc = 5000;

% 非常重要：显式设置参数范围，防止 branch_follow2 自动只扫很短一段
ParamMin = Omega_End;
ParamMax = Omega_Start;

fprintf('\n========== 延拓设置 ==========\n');
fprintf('Omega_Start = %.6f\n', Omega_Start);
fprintf('Omega_End   = %.6f\n', Omega_End);
fprintf('Omega_Step  = %.6f\n', Omega_Step);
fprintf('ParamMin    = %.6f\n', ParamMin);
fprintf('ParamMax    = %.6f\n', ParamMax);

%% -----------------------------
% 4) 构造前两个初始点
%% -----------------------------
fprintf('\n========== 构造初始点 ==========\n');

y_init0 = zeros(16,1);
y_init0(16) = Omega_Start;

[x0_full, ok0, R0n] = newton('nondim_temp2', y_init0, sysP);

fprintf('Point 0: Omega=%.6f, ok=%d, R=%.3e\n', Omega_Start, ok0, R0n);

if ~ok0 || R0n > 1e-6
    error('第一个高频初始点 Newton 未可靠收敛。');
end

x0 = x0_full(1:15);

Omega_Next = Omega_Start + Omega_Step;

y_init1 = [x0; Omega_Next];

[x1_full, ok1, R1n] = newton('nondim_temp2', y_init1, sysP);

fprintf('Point 1: Omega=%.6f, ok=%d, R=%.3e\n', Omega_Next, ok1, R1n);

if ~ok1 || R1n > 1e-6
    error('第二个初始点 Newton 未可靠收敛。');
end

x1 = x1_full(1:15);

%% -----------------------------
% 5) 调用 branch_follow2
%% -----------------------------
fprintf('\n========== 开始 branch_follow2 ==========\n');

[x_branch, conv_flag] = branch_follow2('nondim_temp2', ...
                                       nStepsArc, ...
                                       Omega_Start, ...
                                       Omega_Next, ...
                                       x0, x1, sysP);

fprintf('\n========== branch_follow2 完成 ==========\n');
fprintf('\n========== branch_follow2 完成 ==========\n');

if isstruct(conv_flag)
    fprintf('conv_flag 是 struct，内容如下：\n');
    disp(conv_flag);
elseif ischar(conv_flag) || isstring(conv_flag)
    fprintf('conv_flag = %s\n', string(conv_flag));
elseif isnumeric(conv_flag) || islogical(conv_flag)
    fprintf('conv_flag = %g\n', conv_flag);
else
    fprintf('conv_flag 类型未知：%s\n', class(conv_flag));
    disp(conv_flag);
end

if size(x_branch,1) ~= 16
    error('branch_follow2 输出维度不是 16 x N。');
end

%% -----------------------------
% 6) 检查所有延拓点残差
%% -----------------------------
fprintf('\n========== 全分支残差检查 ==========\n');

Npts = size(x_branch,2);
Rnorm_all = zeros(1,Npts);

for k = 1:Npts
    yk = x_branch(:,k);
    Rk = nondim_temp2(yk, sysP);
    Rnorm_all(k) = norm(Rk);
end

fprintf('Residual max    = %.6e\n', max(Rnorm_all));
fprintf('Residual median = %.6e\n', median(Rnorm_all));
fprintf('Residual mean   = %.6e\n', mean(Rnorm_all));

bad_idx = find(Rnorm_all > 1e-5);

if isempty(bad_idx)
    fprintf('判定：所有延拓点残差均小于 1e-5，通过。\n');
else
    fprintf('警告：有 %d 个点残差 > 1e-5。\n', numel(bad_idx));
    fprintf('最大残差点 index=%d, Omega=%.6f, R=%.3e\n', ...
            bad_idx(1), x_branch(16,bad_idx(1)), Rnorm_all(bad_idx(1)));
end

%% -----------------------------
% 7) 抽样点 fixed-Omega Newton 复核
%% -----------------------------
fprintf('\n========== 抽样点 fixed-Omega Newton 复核 ==========\n');

% 抽样数量
nSample = 10;

% 避免取首尾过于靠近边界的点
if Npts <= nSample
    sample_idx = 1:Npts;
else
    sample_idx = unique(round(linspace(1, Npts, nSample)));
end

sample_results = struct([]);

for ii = 1:numel(sample_idx)

    idx = sample_idx(ii);

    y_branch = x_branch(:,idx);
    Omega_i  = y_branch(16);

    % 用 branch 点作为初值，固定 Omega_i 重新 Newton
    [y_fix, ok_fix, R_fix] = newton('nondim_temp2', y_branch, sysP);

    % 状态误差
    state_err = norm(y_fix(1:15) - y_branch(1:15)) / ...
                max(1e-14, norm(y_fix(1:15)));

    % 传递力幅值和 dB
    TF_branch = calc_TF_point(y_branch, sysP, Fw);
    TF_fix    = calc_TF_point(y_fix,    sysP, Fw);

    TF_abs_err = abs(TF_branch - TF_fix);
    TF_rel_err = TF_abs_err / max(1e-14, abs(TF_fix));

    TFdB_branch = 20*log10(max(TF_branch, 1e-300));
    TFdB_fix    = 20*log10(max(TF_fix,    1e-300));
    TFdB_err    = abs(TFdB_branch - TFdB_fix);

    sample_results(ii).idx = idx;
    sample_results(ii).Omega = Omega_i;
    sample_results(ii).ok = ok_fix;
    sample_results(ii).R = R_fix;
    sample_results(ii).state_err = state_err;
    sample_results(ii).TF_rel_err = TF_rel_err;
    sample_results(ii).TFdB_err = TFdB_err;

    fprintf(['idx=%5d | Omega=%9.5f | ok=%d | R=%9.2e | ', ...
             'state_err=%9.2e | TF_rel=%9.2e | TFdB_err=%9.2e dB\n'], ...
             idx, Omega_i, ok_fix, R_fix, ...
             state_err, TF_rel_err, TFdB_err);
end

%% -----------------------------
% 8) 汇总判据
%% -----------------------------
R_sample      = [sample_results.R];
state_err_all = [sample_results.state_err];
TFrel_all     = [sample_results.TF_rel_err];
TFdBerr_all   = [sample_results.TFdB_err];

fprintf('\n========== 抽样复核汇总 ==========\n');
fprintf('sample residual max  = %.6e\n', max(R_sample));
fprintf('state_err max        = %.6e\n', max(state_err_all));
fprintf('TF_rel_err max       = %.6e\n', max(TFrel_all));
fprintf('TFdB_err max         = %.6e dB\n', max(TFdBerr_all));

if max(R_sample) < 1e-6 && max(state_err_all) < 1e-5 && max(TFdBerr_all) < 1e-4
    fprintf('判定：branch_follow2 抽样复核很好，分支点与定频 Newton 高度一致。\n');
elseif max(R_sample) < 1e-5 && max(state_err_all) < 1e-4 && max(TFdBerr_all) < 1e-3
    fprintf('判定：branch_follow2 抽样复核通过，误差在可接受范围内。\n');
else
    fprintf('判定：branch_follow2 需要进一步检查，可能存在飞点、残差偏大或重复分支问题。\n');
end

%% -----------------------------
% 9) 绘图：残差与 TF
%% -----------------------------
Omega_all = x_branch(16,:);

TF_all = zeros(1,Npts);
for k = 1:Npts
    TF_all(k) = calc_TF_point(x_branch(:,k), sysP, Fw);
end
TFdB_all = 20*log10(max(TF_all, 1e-300));

figure('Color','w','Position',[100 100 1200 450]);

subplot(1,2,1);
semilogy(Omega_all, Rnorm_all, 'b.-', 'LineWidth', 1.2);
grid on;
xlabel('\Omega');
ylabel('||R||');
title('Residual along branch');

subplot(1,2,2);
plot(Omega_all, TFdB_all, 'k-', 'LineWidth', 1.5); hold on;
plot(Omega_all(sample_idx), TFdB_all(sample_idx), 'ro', 'MarkerSize', 7);
grid on;
xlabel('\Omega');
ylabel('TF (dB)');
title('Branch TF and sampled validation points');
legend('Branch','Sampled points','Location','best');

fprintf('\n========== branch_follow2 验证完成 ==========\n');

%% =========================================================
% 局部函数：计算单点基础端力传递率
%% =========================================================
function TF = calc_TF_point(y, sysP, Fw)

    Om  = y(16);
    be2 = sysP(2);
    mu  = sysP(3);
    ze1 = sysP(6);
    ga2 = sysP(11);

    x2 = y(6:10);

    D = zeros(5);
    D(2,3) = Om;
    D(3,2) = -Om;
    D(4,5) = 3*Om;
    D(5,4) = -3*Om;

    x2_dot = D * x2;
    x2_cub = cubic_proj_013_local(x2);

    ft = be2*x2 + ga2*x2_cub + 2*mu*ze1*x2_dot;

    ft1 = hypot(ft(2), ft(3));
    ft3 = hypot(ft(4), ft(5));

    % 这里采用 1+3 次谐波合成幅值
    ft_amp = hypot(ft1, ft3);

    TF = ft_amp / Fw;
end

%% =========================================================
% 局部函数：三次项 AFT 投影
%% =========================================================
function cubic = cubic_proj_013_local(u)

    N = 64;
    t = (0:N-1)'*(2*pi/N);

    dc = ones(N,1);
    c1 = cos(t);
    s1 = sin(t);
    c3 = cos(3*t);
    s3 = sin(3*t);

    T = [dc, c1, s1, c3, s3];

    Inv = [dc, 2*c1, 2*s1, 2*c3, 2*s3]';
    Tinv = (1/N) * Inv;
    Tinv(1,:) = (1/N) * dc';

    cubic = Tinv * ((T*u).^3);
end