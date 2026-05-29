%% run_op_FRF_only.m
% Only run operator-eliminated model (10 coeff + Omega)
% 输出: 单独的 10 参数模型 FRF 曲线

clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'lib'));
global FixedOmega Fw ParamMin ParamMax

FixedOmega = [];     % sweep Omega
Fw = 0.005;

% Omega range bounds for continuation
ParamMin = 0.1;
ParamMax = 10;

%% ---------- 1) Mechanical parameters ----------
mu   = 0.2;     % 质量比 m2/m1
beta = 2.0;     % 下层竖向线性刚度比
K1   = 1.0;     % 上层水平弹簧刚度比
K2   = 0.5;     % 下层水平弹簧刚度比
U    = 2.0;     % 几何非线性尺度参数
L    = 4/9;     % QZS 长度比

% 反推 v 与非线性系数
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

%% ---------- 2) Electrical parameters ----------
P.lam   = 0.18;
P.kap_e = 0.395;
P.kap_c = 0.032;
P.sigma = 0.623;

% 组装系统参数向量 sysP
sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, ...
        P.lam, P.kap_e, P.kap_c, P.sigma, P.ga2];

%% ---------- 3) Continuation settings ----------
Omega_Start = 10.0;
Omega_Step  = -0.01;
Omega_Next  = Omega_Start + Omega_Step;
Nsteps_op   = 1200;

fprintf('Running operator-eliminated FRF...\n');
fprintf('lam=%.4f, kap_e=%.4f, kap_c=%.4f, sigma=%.4f\n', ...
        P.lam, P.kap_e, P.kap_c, P.sigma);

%% ---------- 4) Initial two points ----------
% 第一个点
y0o = [zeros(10,1); Omega_Start];
[x0o, oko, R0o] = newton('nondim_temp2_op', y0o, sysP);
if ~oko
    error('Op Newton fail @Omega_Start, R=%.2e', R0o);
end

% 第二个点
y1o = [x0o(1:10); Omega_Next];
[x1o, ok1o, R1o] = newton('nondim_temp2_op', y1o, sysP);
if ~ok1o
    error('Op Newton fail @Omega_Next, R=%.2e', R1o);
end

%% ---------- 5) Arc-length continuation ----------
[xres_op, ~] = branch_follow2N('nondim_temp2_op', ...
    Nsteps_op, Omega_Start, Omega_Next, x0o, x1o, sysP);

Om_op = xres_op(11,:).';
TF_op_dB = post_TF_op(Om_op, xres_op, sysP, Fw);

%% ---------- 6) Valid points ----------
validO = isfinite(Om_op) & isfinite(TF_op_dB) & (Om_op > 0);
Om_plot = Om_op(validO);
TF_plot = TF_op_dB(validO);

%% ---------- 7) Plot ----------
figure('Color','w','Position',[120 120 860 420]);
hold on; grid on; box on;
set(gca,'XScale','log');

plot(Om_plot, TF_plot, 'b-', 'LineWidth', 1.6);

xlabel('\Omega (log)');
ylabel('TF (dB)');
title('Operator-eliminated model FRF (10 coefficients)');
xlim([ParamMin, ParamMax]);

%% ---------- 8) Optional data export ----------
FRF_op = [Om_plot, TF_plot];
assignin('base', 'FRF_op', FRF_op);
assignin('base', 'Om_op', Om_plot);
assignin('base', 'TF_op_dB', TF_plot);

fprintf('Done. %d valid points obtained.\n', numel(Om_plot));

%% ===================== Helpers =====================
function TF_dB = post_TF_op(Om, xres, sysP, Fw)
    be2 = sysP(2);
    mu  = sysP(3);
    ze2 = sysP(6);
    ga2 = sysP(11);

    % 对于 10 参数模型，下层位移谐波系数仍取第 6:10 行
    % [x20, a21, b21, a23, b23]
    x2 = xres(6:10,:).';

    % 谐波导数
    x2_dot = zeros(size(x2));
    x2_dot(:,2) = Om .* x2(:,3);
    x2_dot(:,3) = -Om .* x2(:,2);
    x2_dot(:,4) = 3*Om .* x2(:,5);
    x2_dot(:,5) = -3*Om .* x2(:,4);

    % 三次项投影
    x2_cub = cubic_proj_013_batch(x2);

    % 传递力
    ft = be2*x2 + ga2*x2_cub + 2*mu*ze2*x2_dot;

    % 基波 + 三次谐波幅值合成
    ft1 = hypot(ft(:,2), ft(:,3));
    ft3 = hypot(ft(:,4), ft(:,5));
    ft_amp = hypot(ft1, ft3);

    TF_dB = 20*log10(max(ft_amp ./ Fw, 1e-300));
end

function cubic = cubic_proj_013_batch(U)
    [~, T_mat, T_inv] = get_AFT_matrices_local();
    cubic = (T_inv * ((T_mat * U.').^3)).';
end

function [N, T_mat, T_inv] = get_AFT_matrices_local()
    persistent pN pT pTinv
    if isempty(pN)
        pN = 64;
        t = (0:pN-1)'*(2*pi/pN);
        c1 = cos(t); s1 = sin(t);
        c3 = cos(3*t); s3 = sin(3*t);
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