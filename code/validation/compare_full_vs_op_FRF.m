%% compare_full_vs_op_FRF.m
% Compare full model (15 coeff) vs operator-eliminated model (10 coeff)

clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'lib'));
global FixedOmega Fw ParamMin ParamMax

FixedOmega = [];          % sweep Omega
Fw = 0.005;

% Omega range bounds for continuation
ParamMin = 0.1;
ParamMax = 10;

mu   = 0.2;     % 质量比 m2/m1
beta = 2.0;     % 下层竖向线性刚度比
K1   = 1.0;     % 上层水平弹簧刚度比
K2   = 0.5;     % 下层水平弹簧刚度比
U    = 2.0;     % 几何非线性尺度参数
L    = 4/9;     % QZS 长度比

% 反推 v 与非线性系数
v = 2.5;        % 由 L=4/9, K1=1 反推
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
P.ze1 = 0.05;   % 下层阻尼比

% 待验证的电路参数（原第二组参数）
P.lam   = 0.18;
P.kap_e = 0.395;
P.kap_c = 0.032;
P.sigma = 0.623;

% 组装系统参数向量 sysP
sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, ...
        P.lam, P.kap_e, P.kap_c, P.sigma, P.ga2];
Omega_Start = 10.0;
Omega_Step  = -0.01;
Omega_Next  = Omega_Start + Omega_Step;
Nsteps_full = 1200;
Nsteps_op   = 1200;

%% ---------- 1) Full model (15 coeff + Omega) ----------
y0 = [zeros(15,1); Omega_Start];
[x0_full, ok0, R0] = newton('nondim_temp2', y0, sysP);
if ~ok0, error('Full Newton fail @Omega_Start, R=%.2e', R0); end

y1 = [x0_full(1:15); Omega_Next];
[x1_full, ok1, R1] = newton('nondim_temp2', y1, sysP);
if ~ok1, error('Full Newton fail @Omega_Next, R=%.2e', R1); end

[xres_full, ~] = branch_follow2('nondim_temp2', Nsteps_full, Omega_Start, Omega_Next, x0_full(1:15), x1_full(1:15), sysP);
Om_full = xres_full(16,:).';
TF_full_dB = post_TF(Om_full, xres_full, sysP, Fw, 'full');

%% ---------- 2) Operator-eliminated model (10 coeff + Omega) ----------
y0o = [zeros(10,1); Omega_Start];
[x0o, oko, R0o] = newton('nondim_temp2_op', y0o, sysP);
if ~oko, error('Op Newton fail @Omega_Start, R=%.2e', R0o); end

y1o = [x0o(1:10); Omega_Next];
[x1o, ok1o, R1o] = newton('nondim_temp2_op', y1o, sysP);
if ~ok1o, error('Op Newton fail @Omega_Next, R=%.2e', R1o); end

[xres_op, ~] = branch_follow2N('nondim_temp2_op', Nsteps_op, Omega_Start, Omega_Next, x0o, x1o, sysP);
Om_op = xres_op(11,:).';
TF_op_dB = post_TF(Om_op, xres_op, sysP, Fw, 'op');

%% ---------- 3) Plot overlay ----------
validF = isfinite(Om_full) & isfinite(TF_full_dB) & Om_full>0;
validO = isfinite(Om_op)   & isfinite(TF_op_dB)   & Om_op>0;

figure('Color','w','Position',[120 120 860 420]); hold on; grid on; box on;
set(gca,'XScale','log');
plot(Om_full(validF), TF_full_dB(validF), 'LineWidth', 1.6);
plot(Om_op(validO),   TF_op_dB(validO),   '--', 'LineWidth', 1.6);
xlabel('\Omega (log)'); ylabel('TF (dB)');
legend('Full model (15)','Operator-eliminated (10)','Location','best');
title('FRF comparison');

%% ---------- 4) Metrics table ----------
M = tf_metrics(Om_full(validF), TF_full_dB(validF), Om_op(validO), TF_op_dB(validO));
disp(M);

% Make a small table
T = table(M.E2, M.Einf, M.dTFpk_dB, M.dOmpk_rel, M.dOm0, M.dOm40, ...
    'VariableNames', {'E2_lin','Einf_lin','dTFpk_dB','dOmpk_rel','dOm0','dOm40'});
disp(T);

%% ---------- helpers ----------
function TF_dB = post_TF(Om, xres, sysP, Fw, mode)
    be2 = sysP(2); mu = sysP(3); ze2 = sysP(6); ga2 = sysP(11);
    if strcmp(mode,'full')
        x2 = xres(6:10,:).';
    else
        x2 = xres(6:10,:).';
    end
    x2_dot = zeros(size(x2));
    x2_dot(:,2) = Om .* x2(:,3); x2_dot(:,3) = -Om .* x2(:,2);
    x2_dot(:,4) = 3*Om .* x2(:,5); x2_dot(:,5) = -3*Om .* x2(:,4);

    x2_cub = cubic_proj_013_batch(x2);
    ft = be2*x2 + ga2*x2_cub + 2*mu*ze2*x2_dot;
    ft_amp = hypot(hypot(ft(:,2), ft(:,3)), hypot(ft(:,4), ft(:,5)));
    TF_dB = 20*log10(max(ft_amp ./ Fw, 1e-300));
end

function cubic = cubic_proj_013_batch(U)
    [~, T_mat, T_inv] = get_AFT_matrices_local();
    cubic = (T_inv * ( (T_mat * U.').^3 )).';
end
function [N, T_mat, T_inv] = get_AFT_matrices_local()
    persistent pN pT pTinv
    if isempty(pN)
        pN = 64; t = (0:pN-1)'*(2*pi/pN);
        c1=cos(t); s1=sin(t); c3=cos(3*t); s3=sin(3*t); dc=ones(pN,1);
        pT = [dc, c1, s1, c3, s3];
        Inv = [dc, 2*c1, 2*s1, 2*c3, 2*s3]';
        pTinv = (1/pN) * Inv; pTinv(1,:) = (1/pN) * dc';
    end
    N = pN; T_mat = pT; T_inv = pTinv;
end

function M = tf_metrics(Om_ref, TFdB_ref, Om_pro, TFdB_pro)
    Om_ref = Om_ref(:); TFdB_ref = TFdB_ref(:);
    Om_pro = Om_pro(:); TFdB_pro = TFdB_pro(:);
    [Om_ref, idxr] = sort(Om_ref); TFdB_ref = TFdB_ref(idxr);
    [Om_pro, idxp] = sort(Om_pro); TFdB_pro = TFdB_pro(idxp);

    TFdB_pro_i = interp1(Om_pro, TFdB_pro, Om_ref, 'pchip', 'extrap');

    TF_ref = 10.^(TFdB_ref/20);
    TF_pro = 10.^(TFdB_pro_i/20);

    e = TF_pro - TF_ref;
    M.E2   = norm(e,2) / max(norm(TF_ref,2), 1e-12);
    M.Einf = max(abs(e));

    [M.TFpk_ref_dB, ir] = max(TFdB_ref);
    [M.TFpk_pro_dB, ip] = max(TFdB_pro_i);
    M.Ompk_ref = Om_ref(ir);
    M.Ompk_pro = Om_ref(ip);
    M.dTFpk_dB = M.TFpk_pro_dB - M.TFpk_ref_dB;
    M.dOmpk_rel = (M.Ompk_pro - M.Ompk_ref) / max(M.Ompk_ref,1e-12);

    M.Om0_ref   = crossing_last(Om_ref, TFdB_ref, 0);
    M.Om0_pro   = crossing_last(Om_ref, TFdB_pro_i, 0);
    M.dOm0      = M.Om0_pro - M.Om0_ref;

    M.Om40_ref  = crossing_last(Om_ref, TFdB_ref, -40);
    M.Om40_pro  = crossing_last(Om_ref, TFdB_pro_i, -40);
    M.dOm40     = M.Om40_pro - M.Om40_ref;
end

function w = crossing_last(Om, y, level)
    s = y - level;
    idx = find(s(1:end-1).*s(2:end) <= 0);
    if isempty(idx)
        if all(y < level), w = min(Om);
        else, w = max(Om);
        end
        return;
    end
    k = idx(end);
    w = Om(k) + (Om(k+1)-Om(k))*(0 - s(k))/(s(k+1)-s(k) + 1e-12);
end