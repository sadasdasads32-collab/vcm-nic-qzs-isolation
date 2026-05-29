%% Run_Stability_Boundary_Map.m
% =========================================================================
% 稳定性边界图：在 (Fw, Omega) 平面上画出稳定/不稳定区域
%
% 对两组参数分别计算：
%   1. 优化后的 EMSD 参数（来自 unified_optimization 结果）
%   2. 纯机械基线（lam=0, 无电路）
%
% 方法：逐频点扫力，Newton+HBM 求解 → Floquet max|μ|
% =========================================================================

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'lib'));

%% ===================== 参数设置 =====================

% --- 优化后的电路参数 ---
sigma_opt = 1.1506;
kap_e_opt = 1.5222;
kap_c_opt = 0.5743;

% --- 机械参数（与 unified_optimization 一致）---
mu   = 0.2;
beta = 2.0;
K1   = 1.0;
K2   = 0.2;
U    = 2.0;
Lg   = 4/9;
v    = 2.5;

alpha1 = v    - 2*K1*(1-Lg)/Lg;
alpha2 = beta - 2*K2*(1-Lg)/Lg;
gamma1 = K1/(U^2 * Lg^3);
gamma2 = K2/(U^2 * Lg^3);

lam_phys = 0.18;
theta    = sqrt(max(lam_phys, 0));
ze1      = 0.05;

be1 = 1.0;
al1 = alpha1 - be1;
be2 = alpha2;

% --- 网格 ---
N_Fw   = 20;
N_Om   = 30;
Fw_vec = logspace(log10(0.001), log10(0.05), N_Fw);
Om_vec = logspace(log10(0.2), log10(6.0), N_Om);

% --- Floquet 设置 ---
Nt_floquet = 600;
tol_stable = 1.002;

% --- 全局变量（HBM 需要）---
global FixedOmega Fw

fprintf('========== 稳定性边界图计算 ==========\n');
fprintf('网格: %d Fw x %d Omega = %d 点\n', N_Fw, N_Om, N_Fw*N_Om);
fprintf('Fw 范围: [%.4f, %.4f]\n', Fw_vec(1), Fw_vec(end));
fprintf('Omega 范围: [%.2f, %.2f]\n\n', Om_vec(1), Om_vec(end));

%% ===================== 情况 1: 优化参数 =====================

fprintf('--- 情况 1: 优化 EMSD 参数 ---\n');

sysP_opt = [be1, be2, mu, al1, gamma1, ze1, lam_phys, ...
            kap_e_opt, kap_c_opt, sigma_opt, gamma2];

[TF_map1, maxMu_map1, ok_map1] = compute_stability_map(...
    sysP_opt, Fw_vec, Om_vec, Nt_floquet);

%% ===================== 情况 2: 纯机械基线 =====================

fprintf('\n--- 情况 2: 纯机械基线 (无电路) ---\n');

sysP_base = [be1, be2, mu, al1, gamma1, ze1, 0.0, ...
             0.0, 0.0, 0.0, gamma2];

[TF_map2, maxMu_map2, ok_map2] = compute_stability_map(...
    sysP_base, Fw_vec, Om_vec, Nt_floquet);

%% ===================== 出图 =====================

fontName = 'Times New Roman';
fsLab = 12; fsTit = 13;

[FF, OO] = meshgrid(Fw_vec, Om_vec);
FF = FF'; OO = OO';  % rows=Fw, cols=Om

figure('Color','w','Position',[50 50 1400 560]);

% --- 图 A: 优化参数 TF ---
subplot(2,3,1);
pcolor(FF, OO, log10(max(TF_map1, 0.01)));
shading flat; colorbar;
set(gca, 'XScale', 'log', 'YScale', 'log');
colormap(gca, jet);
xlabel('F_w', 'FontName', fontName, 'FontSize', fsLab);
ylabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
title(sprintf('log_{10} TF  (EMSD: \\sigma=%.2f, \\kappa_e=%.2f, \\kappa_c=%.2f)', ...
    sigma_opt, kap_e_opt, kap_c_opt), 'FontName', fontName, 'FontSize', fsTit);

% --- 图 B: 优化参数 Floquet ---
subplot(2,3,2);
maxMu_plot1 = maxMu_map1;
maxMu_plot1(~ok_map1) = NaN;
pcolor(FF, OO, maxMu_plot1);
shading flat; colorbar; caxis([0.8, 1.3]);
set(gca, 'XScale', 'log', 'YScale', 'log');
hold on;
% 画稳定边界
contour(FF, OO, maxMu_map1 <= tol_stable & ok_map1, [0.5, 0.5], ...
    'k-', 'LineWidth', 1.5);
contour(FF, OO, maxMu_map1 <= 1.0 & ok_map1, [0.5, 0.5], ...
    'k--', 'LineWidth', 1.0);
hold off;
xlabel('F_w', 'FontName', fontName, 'FontSize', fsLab);
ylabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
title('Floquet max|\mu| (EMSD)', 'FontName', fontName, 'FontSize', fsTit);

% --- 图 C: 优化参数 稳定/不稳定 ---
subplot(2,3,3);
stable_map1 = double(maxMu_map1 <= tol_stable & ok_map1);
stable_map1(~ok_map1) = 0.5;  % 未收敛 = 灰色
pcolor(FF, OO, stable_map1);
shading flat;
colormap(gca, [0.85 0.85 0.85; 0.2 0.6 1.0; 1.0 0.3 0.3]);
caxis([0 2]);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('F_w', 'FontName', fontName, 'FontSize', fsLab);
ylabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
title('Stable (blue) / Unstable (red) / NoConv (gray)', ...
    'FontName', fontName, 'FontSize', fsTit);

% --- 图 D: 基线 TF ---
subplot(2,3,4);
pcolor(FF, OO, log10(max(TF_map2, 0.01)));
shading flat; colorbar;
set(gca, 'XScale', 'log', 'YScale', 'log');
colormap(gca, jet);
xlabel('F_w', 'FontName', fontName, 'FontSize', fsLab);
ylabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
title('log_{10} TF  (Baseline: no circuit)', 'FontName', fontName, 'FontSize', fsTit);

% --- 图 E: 基线 Floquet ---
subplot(2,3,5);
maxMu_plot2 = maxMu_map2;
maxMu_plot2(~ok_map2) = NaN;
pcolor(FF, OO, maxMu_plot2);
shading flat; colorbar; caxis([0.8, 1.3]);
set(gca, 'XScale', 'log', 'YScale', 'log');
hold on;
contour(FF, OO, maxMu_map2 <= tol_stable & ok_map2, [0.5, 0.5], ...
    'k-', 'LineWidth', 1.5);
contour(FF, OO, maxMu_map2 <= 1.0 & ok_map2, [0.5, 0.5], ...
    'k--', 'LineWidth', 1.0);
hold off;
xlabel('F_w', 'FontName', fontName, 'FontSize', fsLab);
ylabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
title('Floquet max|\mu| (Baseline)', 'FontName', fontName, 'FontSize', fsTit);

% --- 图 F: 对比 ---
subplot(2,3,6);
hold on; box on;
% 稳定边界 (EMSD)
contour(FF, OO, maxMu_map1 <= tol_stable & ok_map1, [0.5, 0.5], ...
    'b-', 'LineWidth', 2.0);
% 稳定边界 (Baseline)
contour(FF, OO, maxMu_map2 <= tol_stable & ok_map2, [0.5, 0.5], ...
    'r--', 'LineWidth', 2.0);
set(gca, 'XScale', 'log', 'YScale', 'log');
grid on;
xlabel('F_w', 'FontName', fontName, 'FontSize', fsLab);
ylabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
title('Stability Boundary', 'FontName', fontName, 'FontSize', fsTit);
legend({'EMSD (optimized)', 'Baseline (no circuit)'}, ...
    'Location', 'best', 'FontName', fontName, 'FontSize', 10);

sgtitle(sprintf(['Stability Boundary Map in (F_w, \\Omega) Plane\n' ...
    'EMSD: \\sigma=%.4f, \\kappa_e=%.4f, \\kappa_c=%.4f'], ...
    sigma_opt, kap_e_opt, kap_c_opt), ...
    'FontName', fontName, 'FontSize', 14);

%% ===================== 保存 =====================

out_dir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

% 保存数据
save(fullfile(out_dir, ['stability_map_' timestamp '.mat']), ...
    'Fw_vec', 'Om_vec', ...
    'TF_map1', 'maxMu_map1', 'ok_map1', ...
    'TF_map2', 'maxMu_map2', 'ok_map2', ...
    'sigma_opt', 'kap_e_opt', 'kap_c_opt', 'tol_stable');

% 导出 PDF
exportgraphics(gcf, fullfile(out_dir, ['stability_boundary_map_' timestamp '.pdf']), ...
    'ContentType', 'vector');

fprintf('\n数据已保存到: %s\n', out_dir);
fprintf('图片已导出\n');

%% ===================== 子函数 =====================

function [TF_map, maxMu_map, ok_map] = compute_stability_map(...
    sysP, Fw_vec, Om_vec, Nt_floquet)

    global FixedOmega Fw

    N_Fw = length(Fw_vec);
    N_Om = length(Om_vec);
    TF_map    = nan(N_Fw, N_Om);
    maxMu_map = nan(N_Fw, N_Om);
    ok_map    = false(N_Fw, N_Om);

    for iOm = 1:N_Om
        Om = Om_vec(iOm);
        FixedOmega = Om;

        % 从最低 Fw 开始，零初值
        Fw = Fw_vec(1);
        y_guess = [zeros(15,1); Fw];
        first_ok = false;

        for iFw = 1:N_Fw
            Fw = Fw_vec(iFw);

            % Newton 求解
            try
                y_sol = newton('nondim_temp2', y_guess, sysP);
            catch
                % 不收敛：保留 NaN，重置初值继续
                y_guess = [zeros(15,1); Fw];
                continue;
            end

            xc = y_sol(1:15);
            y_guess = [xc; Fw];  % 延续初值
            ok_map(iFw, iOm) = true;
            if ~first_ok, first_ok = true; end

            % TF
            TF_map(iFw, iOm) = compute_TF_fast(xc, sysP, Om, Fw);

            % Floquet
            maxMu_map(iFw, iOm) = compute_floquet_fast(xc, sysP, Om, Nt_floquet);
        end

        if mod(iOm, 5) == 0 || iOm == 1
            fprintf('  Omega=%5.2f (%2d/%d): %d/%d 频点收敛\n', ...
                Om, iOm, N_Om, nnz(ok_map(:,iOm)), N_Fw);
        end
    end
end
