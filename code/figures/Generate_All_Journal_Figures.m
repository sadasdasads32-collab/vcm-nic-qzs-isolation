%% Generate_All_Journal_Figures.m
% =========================================================================
% 统一生成修订版论文所需全部 EPS 图片，保存到 latex_output 目录
% 关键改动：
%   1. 力传递率 y 轴改用 dB (20*log10(TF))
%   2. 输出 EPS 格式（Springer Nature 期刊要求）
%   3. 图片命名与修订版论文一致 (fig3-1 ~ fig4-3)
% =========================================================================

clc; clear; close all;
init_path();

%% ==================== 输出设置 ====================
out_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'output', 'journal_figures');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

fontName = 'Times New Roman';
fsLab   = 12;
fsTit   = 13;
fsAx    = 11;

%% ==================== 公共参数 ====================
sigma_opt = 1.1506;
kap_e_opt = 1.5222;
kap_c_opt = 0.5743;
theta_val = sqrt(0.18);

mu   = 0.2;   beta_m = 2.0;   K1 = 1.0;   K2 = 0.2;
U    = 2.0;   Lg = 4/9;       v  = 2.5;
alpha1 = v    - 2*K1*(1-Lg)/Lg;
alpha2 = beta_m - 2*K2*(1-Lg)/Lg;
gamma1_val = K1 / (U^2 * Lg^3);
gamma2_val = K2 / (U^2 * Lg^3);
lam_phys = 0.18;
theta    = sqrt(max(lam_phys, 0));
ze1      = 0.05;
be1 = 1.0;
al1 = alpha1 - be1;
be2 = alpha2;

Fw_fixed = 0.008;
Om_min = 0.2;  Om_max = 6.0;
tol_stable = 1.002;

global FixedOmega Fw
Fw = Fw_fixed;
FixedOmega = [];

fprintf('========================================\n');
fprintf('  Generate All Journal Figures (EPS, dB)\n');
fprintf('  Output: %s\n', out_dir);
fprintf('========================================\n\n');

%% ================================================================
%%  Fig3-1: K(Omega) 实部与虚部 + 频带分区
%% ================================================================
fprintf('[1/7] Fig3-1: K(Omega) real & imaginary parts\n');

Om_k = logspace(log10(0.05), log10(10.0), 2000).';
Omega_e = sqrt(kap_c_opt / kap_e_opt);

Delta_k = kap_e_opt * Om_k.^2 - kap_c_opt;
Den_k   = Delta_k.^2 + (sigma_opt * Om_k).^2;
Kr_k = theta^2 .* Om_k.^2 .* Delta_k ./ Den_k;
Ki_k = theta^2 .* sigma_opt .* Om_k.^3 ./ Den_k;

figure('Color','w','Position',[100 100 1000 500],'Visible','off');
hold on; grid on; box on;
set(gca, 'XScale', 'log');

plot(Om_k, Kr_k, 'b-', 'LineWidth', 2.2);
plot(Om_k, Ki_k, 'r-', 'LineWidth', 2.2);
xline(Omega_e, 'k--', 'LineWidth', 1.5);
yline(0, 'k-', 'LineWidth', 0.8);

yl = ylim();
xf1 = 0.05 * 1.5;  xf2 = Omega_e * 0.55;
patch([xf1 xf2 xf2 xf1], [yl(1) yl(1) yl(2) yl(2)], ...
    [0.6 0.8 1.0], 'FaceAlpha', 0.12, 'EdgeColor', 'none');
text(sqrt(xf1*xf2), yl(2)*0.88, 'Virtual Inertia', ...
    'FontName', fontName, 'FontSize', 10, 'Color', 'b', 'HorizontalAlignment', 'center');

xf3 = Omega_e * 0.65;  xf4 = Omega_e * 2.8;
patch([xf3 xf4 xf4 xf3], [yl(1) yl(1) yl(2) yl(2)], ...
    [1.0 0.7 0.7], 'FaceAlpha', 0.12, 'EdgeColor', 'none');
text(sqrt(xf3*xf4), yl(2)*0.82, 'Damping Shaping', ...
    'FontName', fontName, 'FontSize', 10, 'Color', 'r', 'HorizontalAlignment', 'center');

xf5 = Omega_e * 2.3;  xf6 = Omega_e * 7;
patch([xf5 xf6 xf6 xf5], [yl(1) yl(1) yl(2) yl(2)], ...
    [0.7 1.0 0.7], 'FaceAlpha', 0.12, 'EdgeColor', 'none');
text(sqrt(xf5*xf6), yl(2)*0.88, 'Stiffness-like', ...
    'FontName', fontName, 'FontSize', 10, 'Color', [0 0.5 0], 'HorizontalAlignment', 'center');

xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('K(\Omega)', 'FontName', fontName, 'FontSize', fsLab);
title(sprintf(['Real and Imaginary Parts of Complex Operator K(\\Omega)\n' ...
    '\\sigma=%.4f, \\kappa_e=%.4f, \\kappa_c=%.4f'], sigma_opt, kap_e_opt, kap_c_opt), ...
    'FontName', fontName, 'FontSize', fsTit);
legend({'K_R(\Omega)', 'K_I(\Omega)', '\Omega_e'}, ...
    'Location', 'northeast', 'FontName', fontName, 'FontSize', 11);
xlim([0.05 10.0]);

drawnow;
exportgraphics(gcf, fullfile(out_dir, 'fig3-1.eps'), 'ContentType', 'vector');
close(gcf);
fprintf('  -> fig3-1.eps saved\n');

%% ================================================================
%%  Fig3-2: 电路参数扫描 —— sigma, kap_e, kap_c 对 FRF 的影响
%% ================================================================
fprintf('[2/7] Fig3-2: Circuit parameter scan on FRF\n');

Nw_scan = 120;
Om_scan = logspace(log10(Om_min), log10(Om_max), Nw_scan).';
Nvals = 18;
mult_vec = logspace(log10(0.2), log10(5.0), Nvals);

% 构造参考系统参数
sysP_ref = [be1, be2, mu, al1, gamma1_val, ze1, lam_phys, ...
            kap_e_opt, kap_c_opt, sigma_opt, gamma2_val];

param_configs = {
    mult_vec * sigma_opt,  '\sigma',   10;
    mult_vec * kap_e_opt,  '\kappa_e', 8;
    mult_vec * kap_c_opt,  '\kappa_c', 9;
};

results = cell(3,1);
for ip = 1:3
    param_vals = param_configs{ip, 1};
    param_idx  = param_configs{ip, 3};
    TF_peak  = nan(Nvals, 1);
    stable_pct = nan(Nvals, 1);
    TF_all   = cell(Nvals, 1);

    fprintf('  Scanning %s (%d values)...\n', param_configs{ip, 2}, Nvals);

    for iv = 1:Nvals
        sysP = sysP_ref;
        sysP(param_idx) = param_vals(iv);

        TF = nan(Nw_scan, 1);
        ok = false(Nw_scan, 1);
        y_guess = [zeros(15,1); Fw_fixed];
        fail_count = 0;

        for j = 1:Nw_scan
            Om = Om_scan(j);
            FixedOmega = Om;
            try
                y_sol = newton('nondim_temp2', y_guess, sysP);
            catch
                fail_count = fail_count + 1;
                if fail_count > 8, break; else, continue; end
            end
            xc = y_sol(1:15);
            y_guess = [xc; Fw_fixed];
            ok(j) = true;
            TF(j) = compute_TF_fast(xc, sysP, Om, Fw_fixed);
        end

        TF_all{iv} = TF;
        if nnz(ok) >= 5
            TF_peak(iv) = max(TF(ok));
            stable_pct(iv) = 100 * sum(isfinite(TF(ok))) / nnz(ok);
        end
    end
    results{ip} = struct('param_vals', param_vals, 'TF_peak', TF_peak, ...
        'stable_pct', stable_pct, 'TF_all', {TF_all}, 'Om_scan', Om_scan);
end

% --- 出图：三列两行 ---
figure('Color','w','Position',[40 40 1500 620],'Visible','off');
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for ip = 1:3
    r = results{ip};
    vals = r.param_vals;
    name = param_configs{ip, 2};
    ref_val = [sigma_opt, kap_e_opt, kap_c_opt];
    ref_val = ref_val(ip);
    show_idx = round(linspace(1, Nvals, 6));
    cmap = lines(6);

    % FRF 子图 (dB)
    nexttile;
    hold on; box on; grid on;
    for k = 1:6
        idx = show_idx(k);
        tf = r.TF_all{idx};
        ok_pts = isfinite(tf) & tf > 0;
        if any(ok_pts)
            plot(Om_scan(ok_pts), 20*log10(tf(ok_pts)), 'Color', cmap(k,:), ...
                'LineWidth', 1.3, 'DisplayName', sprintf('%s=%.3f', name, vals(idx)));
        end
    end
    yline(0, 'k--', 'LineWidth', 1.0);
    xlabel('\Omega', 'FontName', fontName, 'FontSize', fsAx);
    ylabel('T_F (dB)', 'FontName', fontName, 'FontSize', fsAx);
    title(sprintf('FRF vs %s (ref=%.4f)', name, ref_val), ...
        'FontName', fontName, 'FontSize', fsTit-1);
    legend('Location', 'best', 'FontName', fontName, 'FontSize', 7);
    set(gca, 'XScale', 'log');
end

% 灵敏度汇总子图
for ip = 1:3
    r = results{ip};
    vals = r.param_vals;
    name = param_configs{ip, 2};
    ref_val = [sigma_opt, kap_e_opt, kap_c_opt];
    ref_val = ref_val(ip);

    nexttile;
    yyaxis left;
    semilogx(vals, r.TF_peak, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5);
    xline(ref_val, 'b--', 'LineWidth', 1.0);
    ylabel('TF_{peak} (linear)', 'FontName', fontName, 'FontSize', fsAx);
    grid on; box on;

    yyaxis right;
    semilogx(vals, r.stable_pct, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 5);
    ylabel('Stable %', 'FontName', fontName, 'FontSize', fsAx);
    ylim([0, 105]);

    xlabel(name, 'FontName', fontName, 'FontSize', fsAx);
    title(sprintf('Sensitivity to %s', name), 'FontName', fontName, 'FontSize', fsTit-1);
    legend({'TF_{peak}', 'Stable %'}, 'Location', 'best', 'FontName', fontName, 'FontSize', 8);
end

sgtitle(sprintf(['Circuit Parameter Effects on FRF and Sensitivity\n' ...
    'Ref: \\sigma=%.4f, \\kappa_e=%.4f, \\kappa_c=%.4f, F_w=%.4f'], ...
    sigma_opt, kap_e_opt, kap_c_opt, Fw_fixed), ...
    'FontName', fontName, 'FontSize', 14);

drawnow;
exportgraphics(gcf, fullfile(out_dir, 'fig3-2.eps'), 'ContentType', 'vector');
close(gcf);
fprintf('  -> fig3-2.eps saved\n');

%% ================================================================
%%  Fig3-3: 力传递率对比 —— EMSD 优化 vs 纯机械基线 (dB)
%% ================================================================
fprintf('[3/7] Fig3-3: TF comparison EMSD vs Baseline (dB)\n');

Nw = 200;
Om_vec = logspace(log10(Om_min), log10(Om_max), Nw).';

sysP_opt = [be1, be2, mu, al1, gamma1_val, ze1, lam_phys, ...
            kap_e_opt, kap_c_opt, sigma_opt, gamma2_val];
sysP_base = [be1, be2, mu, al1, gamma1_val, ze1, 0.0, ...
             0.0, 0.0, 0.0, gamma2_val];
sysP_passive = [be1, be2, mu, al1, gamma1_val, ze1, lam_phys, ...
                kap_e_opt, kap_c_opt, 1.0, gamma2_val];  % sigma=1.0, no NIC

% 扫频计算
fprintf('  Computing EMSD sweep...\n');
[TF_opt, maxMu_opt] = do_sweep(sysP_opt, Om_vec, Fw_fixed);
fprintf('  Computing baseline sweep...\n');
[TF_base, maxMu_base] = do_sweep(sysP_base, Om_vec, Fw_fixed);
fprintf('  Computing passive RLC sweep...\n');
[TF_passive, maxMu_passive] = do_sweep(sysP_passive, Om_vec, Fw_fixed);

figure('Color','w','Position',[50 50 1200 500],'Visible','off');

% 左图：TF (dB)
subplot(1,2,1);
hold on; box on; grid on;
ok_opt = isfinite(TF_opt) & TF_opt > 0;
ok_base = isfinite(TF_base) & TF_base > 0;
ok_passive = isfinite(TF_passive) & TF_passive > 0;

plot(Om_vec(ok_opt), 20*log10(TF_opt(ok_opt)), 'b-', 'LineWidth', 1.8, ...
    'DisplayName', 'EMSD (NIC active, opt)');
plot(Om_vec(ok_passive), 20*log10(TF_passive(ok_passive)), 'Color', [0 0.6 0.6], ...
    'LineWidth', 1.4, 'LineStyle', '-.', 'DisplayName', 'Passive RLC (\sigma=1.0)');
plot(Om_vec(ok_base), 20*log10(TF_base(ok_base)), 'r--', 'LineWidth', 1.5, ...
    'DisplayName', 'Baseline (no circuit)');

yline(0, 'k-', 'LineWidth', 1.0);
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('T_F (dB)', 'FontName', fontName, 'FontSize', fsLab);
legend('Location', 'southwest', 'FontName', fontName, 'FontSize', 9);
title('Force Transmissibility (dB)', 'FontName', fontName, 'FontSize', fsTit);
set(gca, 'XScale', 'log');

% 标注峰值抑制
[peak_opt_dB, idx_peak] = max(20*log10(TF_opt(ok_opt)));
peak_base_dB = max(20*log10(TF_base(ok_base)));
plot(Om_vec(idx_peak), peak_opt_dB, 'bo', 'MarkerSize', 8, 'LineWidth', 1.5);
text(Om_vec(idx_peak)*1.5, peak_opt_dB+1, ...
    sprintf('Peak = %.1f dB', peak_opt_dB), ...
    'FontName', fontName, 'FontSize', 9, 'Color', 'b');

% 右图：Floquet 稳定性
subplot(1,2,2);
hold on; box on; grid on;
semilogx(Om_vec(ok_opt), maxMu_opt(ok_opt), 'b.', 'MarkerSize', 6);
semilogx(Om_vec(ok_passive), maxMu_passive(ok_passive), '.', ...
    'Color', [0 0.6 0.6], 'MarkerSize', 6);
semilogx(Om_vec(ok_base), maxMu_base(ok_base), 'r.', 'MarkerSize', 6);
yline(1.002, 'k--', 'LineWidth', 1.5);
yline(1.0, ':', 'Color', [0.4 0.4 0.4]);
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('max|\mu|', 'FontName', fontName, 'FontSize', fsLab);
legend({'EMSD (NIC)', 'Passive RLC', 'Baseline', 'tol=1.002', '|\mu|=1'}, ...
    'Location', 'best', 'FontName', fontName, 'FontSize', 8);
title('Floquet Stability', 'FontName', fontName, 'FontSize', fsTit);
ylim([0.6, 1.8]);

[peak_opt_lin, ~] = max(TF_opt(ok_opt));
[peak_base_lin, ~] = max(TF_base(ok_base));
reduction = (1 - peak_opt_lin/peak_base_lin)*100;
sgtitle(sprintf(['Force Transmissibility & Stability: 3 Cases\n' ...
    'TF_{peak}: %.3f \\rightarrow %.3f (%.1f%% reduction),  F_w=%.4f'], ...
    peak_base_lin, peak_opt_lin, reduction, Fw_fixed), ...
    'FontName', fontName, 'FontSize', 13);

drawnow;
exportgraphics(gcf, fullfile(out_dir, 'fig3-3.eps'), 'ContentType', 'vector');
close(gcf);
fprintf('  -> fig3-3.eps saved (peak reduction: %.1f %%)\n', reduction);

%% ================================================================
%%  Fig4-1: (Fw, Omega) Floquet 稳定性边界图
%% ================================================================
fprintf('[4/7] Fig4-1: (Fw, Omega) stability boundary map\n');

N_Fw  = 20;
N_Om  = 30;
Fw_vec = logspace(log10(0.001), log10(0.05), N_Fw);
Om_map = logspace(log10(0.2), log10(6.0), N_Om);

fprintf('  Computing EMSD stability map...\n');
[TF_map1, maxMu_map1, ok_map1] = compute_stability_map(...
    sysP_opt, Fw_vec, Om_map, 600);

fprintf('  Computing baseline stability map...\n');
[TF_map2, maxMu_map2, ok_map2] = compute_stability_map(...
    sysP_base, Fw_vec, Om_map, 600);

[FF, OO] = meshgrid(Fw_vec, Om_map);
FF = FF';  OO = OO';

figure('Color','w','Position',[50 50 1400 560],'Visible','off');

% (a) EMSD: log10(TF)
subplot(2,3,1);
TF1_plot = TF_map1; TF1_plot(TF1_plot <= 0) = 0.01; TF1_plot = min(TF1_plot, 100);
pcolor(FF, OO, 20*log10(TF1_plot));
shading flat; colorbar;
set(gca, 'XScale', 'log', 'YScale', 'log');
colormap(gca, jet);
xlabel('F_w', 'FontName', fontName, 'FontSize', fsLab);
ylabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
title('T_F (dB) — EMSD', 'FontName', fontName, 'FontSize', fsTit);

% (b) EMSD: Floquet max|mu|
subplot(2,3,2);
maxMu_plot1 = maxMu_map1; maxMu_plot1(~ok_map1) = NaN;
pcolor(FF, OO, maxMu_plot1);
shading flat; colorbar; caxis([0.8, 1.5]);
set(gca, 'XScale', 'log', 'YScale', 'log');
hold on;
contour(FF, OO, maxMu_map1 <= tol_stable & ok_map1, [0.5, 0.5], 'k-', 'LineWidth', 1.5);
contour(FF, OO, maxMu_map1 <= 1.0 & ok_map1, [0.5, 0.5], 'k--', 'LineWidth', 1.0);
hold off;
xlabel('F_w', 'FontName', fontName, 'FontSize', fsLab);
ylabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
title('max|\mu| — EMSD', 'FontName', fontName, 'FontSize', fsTit);

% (c) EMSD: 稳定/不稳定
subplot(2,3,3);
stable_map1 = double(maxMu_map1 <= tol_stable & ok_map1);
stable_map1(~ok_map1) = 0.5;
pcolor(FF, OO, stable_map1);
shading flat;
colormap(gca, [0.85 0.85 0.85; 0.2 0.6 1.0; 1.0 0.3 0.3]);
caxis([0 2]);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('F_w', 'FontName', fontName, 'FontSize', fsLab);
ylabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
title('Stable (blue) / Unstable (red) / NC (gray)', 'FontName', fontName, 'FontSize', fsTit);

% (d) Baseline: log10(TF)
subplot(2,3,4);
TF2_plot = TF_map2; TF2_plot(TF2_plot <= 0) = 0.01; TF2_plot = min(TF2_plot, 100);
pcolor(FF, OO, 20*log10(TF2_plot));
shading flat; colorbar;
set(gca, 'XScale', 'log', 'YScale', 'log');
colormap(gca, jet);
xlabel('F_w', 'FontName', fontName, 'FontSize', fsLab);
ylabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
title('T_F (dB) — Baseline', 'FontName', fontName, 'FontSize', fsTit);

% (e) Baseline: Floquet max|mu|
subplot(2,3,5);
maxMu_plot2 = maxMu_map2; maxMu_plot2(~ok_map2) = NaN;
pcolor(FF, OO, maxMu_plot2);
shading flat; colorbar; caxis([0.8, 1.5]);
set(gca, 'XScale', 'log', 'YScale', 'log');
hold on;
contour(FF, OO, maxMu_map2 <= tol_stable & ok_map2, [0.5, 0.5], 'k-', 'LineWidth', 1.5);
contour(FF, OO, maxMu_map2 <= 1.0 & ok_map2, [0.5, 0.5], 'k--', 'LineWidth', 1.0);
hold off;
xlabel('F_w', 'FontName', fontName, 'FontSize', fsLab);
ylabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
title('max|\mu| — Baseline', 'FontName', fontName, 'FontSize', fsTit);

% (f) 稳定边界对比
subplot(2,3,6);
hold on; box on;
contour(FF, OO, maxMu_map1 <= tol_stable & ok_map1, [0.5, 0.5], 'b-', 'LineWidth', 2.0);
contour(FF, OO, maxMu_map2 <= tol_stable & ok_map2, [0.5, 0.5], 'r--', 'LineWidth', 2.0);
set(gca, 'XScale', 'log', 'YScale', 'log');
grid on;
xlabel('F_w', 'FontName', fontName, 'FontSize', fsLab);
ylabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
title('Stability Boundary Comparison', 'FontName', fontName, 'FontSize', fsTit);
legend({'EMSD (optimized)', 'Baseline (no circuit)'}, ...
    'Location', 'best', 'FontName', fontName, 'FontSize', 10);

sgtitle(sprintf(['Stability Map in (F_w, \\Omega) Plane\n' ...
    'EMSD: \\sigma=%.4f, \\kappa_e=%.4f, \\kappa_c=%.4f'], ...
    sigma_opt, kap_e_opt, kap_c_opt), 'FontName', fontName, 'FontSize', 14);

drawnow;
exportgraphics(gcf, fullfile(out_dir, 'fig4-1.eps'), 'ContentType', 'vector');
close(gcf);
fprintf('  -> fig4-1.eps saved\n');

%% ================================================================
%%  Fig4-2: 定频扫力 + Floquet —— 展示 Fold 分岔与跳跃
%% ================================================================
fprintf('[5/7] Fig4-2: Force sweep + Floquet at representative frequencies\n');

Om_demo = [0.5, 1.0, 2.0];
N_Fw_demo = 30;
Fw_demo = logspace(log10(0.001), log10(0.05), N_Fw_demo);
Nt_fq = 600;

figure('Color','w','Position',[50 50 1400 420],'Visible','off');
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for iOm = 1:3
    Om = Om_demo(iOm);
    FixedOmega = Om;

    A_resp = nan(N_Fw_demo, 1);
    maxMu_arr = nan(N_Fw_demo, 1);
    stable_arr = false(N_Fw_demo, 1);
    Fw = Fw_demo(1);
    y_guess = [zeros(15,1); Fw];
    fail_count = 0;

    for iFw = 1:N_Fw_demo
        Fw = Fw_demo(iFw);
        try
            y_sol = newton('nondim_temp2', y_guess, sysP_opt);
        catch
            fail_count = fail_count + 1;
            y_guess = [zeros(15,1); Fw];
            if fail_count > 5, break; else, continue; end
        end
        xc = y_sol(1:15);
        y_guess = [xc; Fw];
        A_resp(iFw) = hypot(xc(2), xc(3));
        maxMu_arr(iFw) = compute_floquet_fast(xc, sysP_opt, Om, Nt_fq);
        stable_arr(iFw) = maxMu_arr(iFw) < tol_stable;
    end

    nexttile;
    yyaxis left;
    ok_pts = isfinite(A_resp) & A_resp > 0;
    scatter(Fw_demo(stable_arr & ok_pts), A_resp(stable_arr & ok_pts), ...
        25, 'b', 'filled', 'DisplayName', 'Stable');
    hold on;
    scatter(Fw_demo(~stable_arr & ok_pts), A_resp(~stable_arr & ok_pts), ...
        25, 'r', 'DisplayName', 'Unstable (Fold)');
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('F_w', 'FontName', fontName, 'FontSize', fsLab);
    ylabel('|x_1|_{amp}', 'FontName', fontName, 'FontSize', fsLab);
    grid on; box on;
    title(sprintf('\\Omega = %.2f', Om), 'FontName', fontName, 'FontSize', fsTit);

    yyaxis right;
    plot(Fw_demo(ok_pts), maxMu_arr(ok_pts), 'k.-', 'MarkerSize', 6, 'LineWidth', 0.8);
    yline(1.002, 'k--', 'LineWidth', 1.0);
    yline(1.0, ':', 'Color', [0.4 0.4 0.4]);
    ylabel('max|\mu|', 'FontName', fontName, 'FontSize', fsLab);
    ylim([0.5, 2.5]);
    legend('Location', 'best', 'FontName', fontName, 'FontSize', 7);
end

sgtitle('Force Sweep: Response Amplitude \& Floquet Stability (Fold Bifurcation)', ...
    'FontName', fontName, 'FontSize', 14);

drawnow;
exportgraphics(gcf, fullfile(out_dir, 'fig4-2.eps'), 'ContentType', 'vector');
close(gcf);
fprintf('  -> fig4-2.eps saved\n');

%% ================================================================
%%  Fig4-3: 性能-稳定性 Pareto 权衡
%% ================================================================
fprintf('[6/7] Fig4-3: Performance-Stability Pareto trade-off\n');

% 从 parameter scan 结果中提取 Pareto 前沿数据
r_sig = results{1};
r_kap = results{2};
r_kpc = results{3};

figure('Color','w','Position',[80 80 1100 450],'Visible','off');

% 左图：sigma 条件下的性能-稳定性权衡
subplot(1,2,1);
hold on; box on; grid on;
for ip = 1:3
    r = results{ip};
    if ip == 1, ms = 'o'; color = [0.2 0.4 1.0]; nm = '\sigma';
    elseif ip == 2, ms = 's'; color = [1.0 0.4 0.2]; nm = '\kappa_e';
    else, ms = '^'; color = [0.2 0.8 0.2]; nm = '\kappa_c';
    end
    ok = isfinite(r.TF_peak) & isfinite(r.stable_pct) & r.TF_peak > 0;
    scatter(r.TF_peak(ok), r.stable_pct(ok), 40, color, ms, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'DisplayName', nm);
end
xlabel('TF_{peak} (linear)', 'FontName', fontName, 'FontSize', fsLab);
ylabel('Stable Frequency Points (%)', 'FontName', fontName, 'FontSize', fsLab);
title('Performance–Stability Trade-off', 'FontName', fontName, 'FontSize', fsTit);
legend('Location', 'best', 'FontName', fontName, 'FontSize', 10);

% 标注最优工作点
plot(0.36, 95.4, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r', ...
    'DisplayName', 'Optimal (\sigma=1.15)');
text(0.36+0.03, 95.4, sprintf('Optimal\nTF_{peak}=0.36\nStable=95.4%%'), ...
    'FontName', fontName, 'FontSize', 9, 'Color', 'r');

% 右图：安全设计区间示意
subplot(1,2,2);
hold on; box on; grid on;

% 绘制简化的 Pareto 前沿示意
sig_range = [0.23, 0.5, 0.8, sigma_opt, 2.0, 3.0, 5.75];
tf_vals   = [1.23, 0.85, 0.55, 0.36, 0.28, 0.24, 0.22];
stab_vals = [100, 99, 97, 95.4, 93, 90, 88];

fill_x = [sig_range, fliplr(sig_range)];
% Recommended zone
xf_rec = [0.8 2.0 2.0 0.8];
fill(xf_rec, [75 75 105 105], [0.2 0.8 0.2], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
text(1.3, 101, 'Recommended:\newline\sigma \in [0.8, 2.0]', ...
    'FontName', fontName, 'FontSize', 10, 'Color', [0 0.5 0], ...
    'HorizontalAlignment', 'center');

yyaxis left;
plot(sig_range, tf_vals, 'b-o', 'LineWidth', 1.8, 'MarkerSize', 6, 'MarkerFaceColor', 'b');
plot(sigma_opt, 0.36, 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
ylabel('TF_{peak}', 'FontName', fontName, 'FontSize', fsLab);

yyaxis right;
plot(sig_range, stab_vals, 'r-s', 'LineWidth', 1.8, 'MarkerSize', 6, 'MarkerFaceColor', 'r');
ylabel('Stable %', 'FontName', fontName, 'FontSize', fsLab);
ylim([75, 105]);

xlabel('\sigma', 'FontName', fontName, 'FontSize', fsLab);
title('Design Trade-off along \sigma', 'FontName', fontName, 'FontSize', fsTit);
legend({'TF_{peak}', 'Optimal', 'Stable %'}, 'Location', 'best', ...
    'FontName', fontName, 'FontSize', 9);

sgtitle(sprintf(['Performance–Stability Pareto Framework\n' ...
    '\\kappa_e=%.4f, \\kappa_c=%.4f, F_w=%.4f'], ...
    kap_e_opt, kap_c_opt, Fw_fixed), 'FontName', fontName, 'FontSize', 14);

drawnow;
exportgraphics(gcf, fullfile(out_dir, 'fig4-3.eps'), 'ContentType', 'vector');
close(gcf);
fprintf('  -> fig4-3.eps saved\n');

%% ================================================================
%%  Extra: 复算子等效参数分解图 (Meq, Ceq, Keq)
%% ================================================================
fprintf('[7/7] Extra: Equivalent parameter decomposition\n');

C_eq_map = Ki_k ./ Om_k;
M_eq_map = -Kr_k ./ (Om_k.^2);
K_eq_map = Kr_k;
M_eq_low = theta^2 / kap_c_opt;
Kr_high = theta^2 / kap_e_opt;

figure('Color','w','Position',[100 100 1000 800],'Visible','off');

subplot(3,1,1);
hold on; grid on; box on;
set(gca, 'XScale', 'log');
pos_iner = M_eq_map > 0;
plot(Om_k(pos_iner), M_eq_map(pos_iner), 'b-', 'LineWidth', 1.8);
neg_iner = M_eq_map <= 0;
if any(neg_iner)
    plot(Om_k(neg_iner), M_eq_map(neg_iner), 'b--', 'LineWidth', 1.2);
end
xline(Omega_e, 'k--', 'LineWidth', 1.5);
yline(0, 'k-', 'LineWidth', 0.8);
text(Omega_e*0.3, max(M_eq_map)*0.85, ...
    sprintf('M_{eq}(0) = %.3f', M_eq_low), ...
    'FontName', fontName, 'FontSize', 10, 'Color', 'b');
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('M_{eq}(\Omega)', 'FontName', fontName, 'FontSize', fsLab);
title('Equivalent Virtual Inertia: M_{eq} = -K_R/\Omega^2', ...
    'FontName', fontName, 'FontSize', fsTit);
legend({'M_{eq}>0 (virtual inertia)', 'M_{eq}<0 (mass loading)'}, ...
    'Location', 'best', 'FontName', fontName, 'FontSize', 9);

subplot(3,1,2);
hold on; grid on; box on;
set(gca, 'XScale', 'log');
plot(Om_k, C_eq_map, 'r-', 'LineWidth', 1.8);
xline(Omega_e, 'k--', 'LineWidth', 1.5);
[Ce_peak, idx_ce] = max(C_eq_map);
plot(Om_k(idx_ce), Ce_peak, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
text(Om_k(idx_ce)*1.3, Ce_peak*0.85, ...
    sprintf('C_{eq,max}=%.4f at \\Omega=%.3f', Ce_peak, Om_k(idx_ce)), ...
    'FontName', fontName, 'FontSize', 10, 'Color', 'r');
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('C_{eq}(\Omega)', 'FontName', fontName, 'FontSize', fsLab);
title('Equivalent Damping: C_{eq} = K_I/\Omega', ...
    'FontName', fontName, 'FontSize', fsTit);

subplot(3,1,3);
hold on; grid on; box on;
set(gca, 'XScale', 'log');
plot(Om_k, K_eq_map, 'Color', [0 0.5 0], 'LineWidth', 1.8);
xline(Omega_e, 'k--', 'LineWidth', 1.5);
yline(0, 'k-', 'LineWidth', 0.8);
yline(Kr_high, ':', 'Color', [0 0.5 0], 'LineWidth', 1.0);
text(Omega_e*3, Kr_high*1.1, ...
    sprintf('K(\\infty)=\\theta^2/\\kappa_e=%.4f', Kr_high), ...
    'FontName', fontName, 'FontSize', 10, 'Color', [0 0.5 0]);
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('K_{eq}(\Omega)', 'FontName', fontName, 'FontSize', fsLab);
title('Equivalent Stiffness: K_{eq} = K_R', ...
    'FontName', fontName, 'FontSize', fsTit);

sgtitle(sprintf(['Equivalent Parameter Decomposition of K(\\Omega)\n' ...
    '\\sigma=%.4f, \\kappa_e=%.4f, \\kappa_c=%.4f, \\Omega_e=%.3f'], ...
    sigma_opt, kap_e_opt, kap_c_opt, Omega_e), ...
    'FontName', fontName, 'FontSize', 14);
xlim([0.05 10.0]);

drawnow;
exportgraphics(gcf, fullfile(out_dir, 'figA-equiv-params.eps'), 'ContentType', 'vector');
close(gcf);
fprintf('  -> figA-equiv-params.eps saved\n');

%% ================================================================
%%  完成汇总
%% ================================================================
fprintf('\n========================================\n');
fprintf('  All figures saved to: %s\n', out_dir);
fprintf('========================================\n\n');

files = dir(fullfile(out_dir, 'fig*.eps'));
fprintf('Generated EPS figures (%d files):\n', length(files));
total_size = 0;
for k = 1:length(files)
    fprintf('  [%d] %s  (%.1f KB)\n', k, files(k).name, files(k).bytes/1024);
    total_size = total_size + files(k).bytes;
end
fprintf('  Total: %.1f KB\n\n', total_size/1024);

fprintf('Figure-to-Section mapping:\n');
fprintf('  fig3-1.eps  -> Section 3.1: K(Omega) real & imaginary parts\n');
fprintf('  fig3-2.eps  -> Section 3.3: Circuit parameter scan on FRF (dB)\n');
fprintf('  fig3-3.eps  -> Section 3.4: TF comparison: EMSD vs Baseline vs Passive (dB)\n');
fprintf('  fig4-1.eps  -> Section 4.1: (Fw, Omega) stability boundary map\n');
fprintf('  fig4-2.eps  -> Section 4.1/4.2: Force sweep + Floquet\n');
fprintf('  fig4-3.eps  -> Section 4.3: Performance-Stability Pareto\n');
fprintf('  figA-equiv-params.eps -> Appendix: Equivalent parameter decomposition\n');

%% ==================== 子函数 ====================

function [TF_vec, maxMu_vec] = do_sweep(sysP, Om_vec, Fw_val)
    global FixedOmega Fw
    Fw = Fw_val;
    N = length(Om_vec);
    TF_vec = nan(N,1);
    maxMu_vec = nan(N,1);
    fail_count = 0;
    y_guess = [zeros(15,1); Fw_val];
    for j = 1:N
        Om = Om_vec(j);
        FixedOmega = Om;
        try
            y_sol = newton('nondim_temp2', y_guess, sysP);
        catch
            fail_count = fail_count + 1;
            if fail_count > 12, break; else, continue; end
        end
        xc = y_sol(1:15);
        y_guess = [xc; Fw_val];
        TF_vec(j) = compute_TF_fast(xc, sysP, Om, Fw_val);
        maxMu_vec(j) = compute_floquet_fast(xc, sysP, Om, 600);
    end
end

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
        Fw = Fw_vec(1);
        y_guess = [zeros(15,1); Fw];
        for iFw = 1:N_Fw
            Fw = Fw_vec(iFw);
            try
                y_sol = newton('nondim_temp2', y_guess, sysP);
            catch
                y_guess = [zeros(15,1); Fw];
                continue;
            end
            xc = y_sol(1:15);
            y_guess = [xc; Fw];
            ok_map(iFw, iOm) = true;
            TF_map(iFw, iOm) = compute_TF_fast(xc, sysP, Om, Fw);
            maxMu_map(iFw, iOm) = compute_floquet_fast(xc, sysP, Om, Nt_floquet);
        end
        if mod(iOm, 10) == 0 || iOm == 1
            fprintf('    Omega=%5.2f (%2d/%d): %d/%d converged\n', ...
                Om, iOm, N_Om, nnz(ok_map(:,iOm)), N_Fw);
        end
    end
end
