%% 纯电阻电路多组参数对比（Rt = 3, 6, 9, 12）
% 调用 nondim_temp2_1，其余代码结构基本不变

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'lib'));

global Fw FixedOmega
FixedOmega = [];

%% ------------------------------------------------------------
% 1. 物理参数（固定）
% ------------------------------------------------------------
Kt = 7.474;
Ke = 7.474;
m1 = 2.2;
k1 = 3000;
R0 = 3.8;

wn = sqrt(k1/m1);

%% ------------------------------------------------------------
% 2. 基础结构参数（固定）
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

% 固定机电耦合与零电感电容
P.lam   = Kt*Ke*wn/(k1*R0);
P.kap_e = 0.0;
P.kap_c = 0.0;

%% ------------------------------------------------------------
% 3. 激励和延拓设置
% ------------------------------------------------------------
Fw = 0.005;
Omega_Start = 10.0;
Omega_Step  = -0.01;
Omega_Next  = Omega_Start + Omega_Step;

%% ------------------------------------------------------------
% 4. 要对比的 Rt 值
% ------------------------------------------------------------
Rt_values = [3, 6, 9, 12];
n_Rt = length(Rt_values);

% 预存结果
Om_cell  = cell(n_Rt, 1);
TF_cell  = cell(n_Rt, 1);
legend_str = cell(n_Rt, 1);

%% ------------------------------------------------------------
% 5. 循环计算不同 Rt
% ------------------------------------------------------------
for idx = 1:n_Rt
    Rt = Rt_values(idx);
    P.sigma = Rt / R0;
    
    sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, ...
            P.lam, P.kap_e, P.kap_c, P.sigma, P.ga2];
    
    fprintf('\n========== Rt = %.6g (sigma = %.6g) ==========\n', Rt, P.sigma);
    
    % 初始猜测
    y_init = zeros(15,1);
    y_init(end+1) = Omega_Start;
    
    [x0_full, ok0] = newton('nondim_temp2_1', y_init, sysP);
    if ~ok0
        warning('Rt=%.6g: 高频起点求解失败，跳过。', Rt);
        continue;
    end
    x0 = x0_full(1:15);
    
    y_init2 = [x0; Omega_Next];
    [x1_full, ok1] = newton('nondim_temp2_1', y_init2, sysP);
    if ~ok1
        warning('Rt=%.6g: 第二个初始点求解失败，跳过。', Rt);
        continue;
    end
    x1 = x1_full(1:15);
    
    % 弧长延拓
    n_steps = 3000;
    [x_res, ~] = branch_follow2('nondim_temp2_1', ...
                                n_steps, ...
                                Omega_Start, ...
                                Omega_Next, ...
                                x0, ...
                                x1, ...
                                sysP);
    
    % 计算力传递率
    Om  = x_res(16,:).';
    
    be2 = P.be2;
    mu  = P.mu;
    ze1 = P.ze1;
    ga2 = P.ga2;
    
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
    
    TF    = ft_amp ./ Fw;
    TF_dB = 20*log10(max(TF, 1e-300));
    
    valid = isfinite(Om) & isfinite(TF_dB) & (Om > 0);
    Om_cell{idx}  = Om(valid);
    TF_cell{idx}  = TF_dB(valid);
    legend_str{idx} = sprintf('R_t = %.6g', Rt);
    
    fprintf('Rt=%.6g: 有效点数 %d\n', Rt, nnz(valid));
end

%% ------------------------------------------------------------
% 6. 绘图（四张图合一）
% ------------------------------------------------------------
figure('Color','w', 'Position',[150 150 760 520]);
ax = gca;
hold(ax,'on'); grid(ax,'on'); box(ax,'on');
set(ax,'XScale','log');

xlabel(ax, '\Omega (log scale)');
ylabel(ax, 'Force Transmissibility 20log_{10}(|f_t|/f) (dB)');
title(ax, 'BG Model with Resistor-only Circuit: \lambda=4.23, various R_t');
yline(ax, 0, 'k--', '0 dB');

% 颜色/线型可自定义
colors = lines(n_Rt);
for idx = 1:n_Rt
    if ~isempty(Om_cell{idx})
        plot(ax, Om_cell{idx}, TF_cell{idx}, ...
             'Color', colors(idx,:), 'LineWidth', 1.5, ...
             'DisplayName', legend_str{idx});
    end
end

legend(ax, 'Location', 'best');
xlim(ax, [0.1, Omega_Start]);
hold(ax,'off');

%% ------------------------------------------------------------
% 7. 保存所有结果（可选）
% ------------------------------------------------------------
results.Rt_values = Rt_values;
results.Om_cell   = Om_cell;
results.TF_cell   = TF_cell;
results.P         = P;
save('result_resistor_multi.mat', 'results');
fprintf('\n结果已保存到 result_resistor_multi.mat\n');

%% ============================================================
% 辅助函数：与原始代码完全相同
% ============================================================
function cubic = cubic_proj_013_batch(U)
    [~, T_mat, T_inv] = get_AFT_matrices_local();
    X_time  = (T_mat * U.').';
    X3_time = X_time.^3;
    cubic   = (T_inv * X3_time.').';
end

function [N, T_mat, T_inv] = get_AFT_matrices_local()
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