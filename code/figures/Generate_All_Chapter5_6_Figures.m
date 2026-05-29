%% Generate_All_Chapter5_6_Figures.m
% =========================================================================
% 批量生成第5章和第6章所需全部图片，保存到 e:\项目1\论文图\
% =========================================================================

clc; clear; close all;
init_path();

out_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'output', 'journal_figures');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

fontName = 'Times New Roman';
fsLab = 12; fsTit = 14;
export_pdf = true;

fprintf('========================================\n');
fprintf('  生成第5章和第6章论文图片\n');
fprintf('========================================\n\n');

%% ==================== 公共参数 ====================

sigma_opt = 1.1506;
kap_e_opt = 1.5222;
kap_c_opt = 0.5743;

mu   = 0.2;   beta = 2.0;   K1 = 1.0;   K2 = 0.2;
U    = 2.0;   Lg   = 4/9;   v  = 2.5;
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

Fw_fixed = 0.008;
Om_min = 0.2;  Om_max = 6.0;
Nw = 200;  % 频点数（更高分辨率用于论文图）

global FixedOmega Fw
Fw = Fw_fixed;
FixedOmega = [];

Om_vec = logspace(log10(Om_min), log10(Om_max), Nw).';

%% ================================================================
%%  图1: 复动力学算子 K(Omega) —— 第5.4节
%% ================================================================
fprintf('--- 图1: 复动力学算子 K(Omega) ---\n');

Om_k = logspace(log10(0.1), log10(10), 500);
K_op = (theta^2 * Om_k.^2) ./ (kap_e_opt * Om_k.^2 - 1i * sigma_opt * Om_k - kap_c_opt);
Kr = real(K_op);
Ki = imag(K_op);

figure('Color','w','Position',[50 50 1100 420]);

subplot(1,2,1);
hold on; box on; grid on;
plot(Om_k, Kr, 'b-', 'LineWidth', 1.8);
plot(Om_k, Ki, 'r-', 'LineWidth', 1.8);
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('K(\Omega)', 'FontName', fontName, 'FontSize', fsLab);
legend({'K_R(\Omega)', 'K_I(\Omega)'}, 'Location', 'best', ...
    'FontName', fontName, 'FontSize', 11);
title('Complex Dynamics Operator', 'FontName', fontName, 'FontSize', fsTit);
set(gca, 'XScale', 'log');

% 等效参数分解
Ceq = Kr;  % Re[K] 提供等效阻尼
Keq = -Om_k(:) .* Ki(:);  % -Omega * Im[K] = 等效刚度
Meq = Ki(:) ./ Om_k(:);  % Im[K]/Omega = 等效质量

subplot(1,2,2);
hold on; box on; grid on;
plot(Om_k, Ceq, 'b-', 'LineWidth', 1.3, 'DisplayName', 'C_{eq} (damping)');
plot(Om_k, Meq, 'r-', 'LineWidth', 1.3, 'DisplayName', 'M_{eq} (inertia)');
plot(Om_k, Keq, 'Color', [0 0.5 0], 'LineWidth', 1.3, 'DisplayName', 'K_{eq} (stiffness)');
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('Equivalent Parameters', 'FontName', fontName, 'FontSize', fsLab);
legend('Location', 'best', 'FontName', fontName, 'FontSize', 9);
title('Equivalent Parameter Splitting', 'FontName', fontName, 'FontSize', fsTit);
set(gca, 'XScale', 'log');

sgtitle(sprintf(['K(\\Omega) = \\theta^2\\Omega^2/(\\kappa_e\\Omega^2 - j\\sigma\\Omega - \\kappa_c)\n' ...
    '\\sigma=%.4f, \\kappa_e=%.4f, \\kappa_c=%.4f'], sigma_opt, kap_e_opt, kap_c_opt), ...
    'FontName', fontName, 'FontSize', 15);

if export_pdf
    exportgraphics(gcf, fullfile(out_dir, 'Fig5_K_operator.pdf'), 'ContentType', 'vector');
    fprintf('  -> Fig5_K_operator.pdf\n');
end

%% ================================================================
%%  图2: 力传递率对比：优化 EMSD vs 纯机械基线 —— 第5.3节
%% ================================================================
fprintf('--- 图2: 力传递率对比 ---\n');

sysP_opt = [be1, be2, mu, al1, gamma1, ze1, lam_phys, ...
            kap_e_opt, kap_c_opt, sigma_opt, gamma2];
sysP_base = [be1, be2, mu, al1, gamma1, ze1, 0.0, ...
             0.0, 0.0, 0.0, gamma2];

% 扫频计算（内联）
fprintf('  计算 EMSD 扫频...\n');
[TF_opt, maxMu_opt] = do_sweep(sysP_opt, Om_vec, Fw_fixed);
fprintf('  计算基线扫频...\n');
[TF_base, maxMu_base] = do_sweep(sysP_base, Om_vec, Fw_fixed);

figure('Color','w','Position',[50 50 1200 500]);

subplot(1,2,1);
hold on; box on; grid on;
ok_opt = isfinite(TF_opt) & TF_opt > 0;
ok_base = isfinite(TF_base) & TF_base > 0;
TF_opt_dB = 20*log10(TF_opt(ok_opt));
TF_base_dB = 20*log10(TF_base(ok_base));
semilogx(Om_vec(ok_opt), TF_opt_dB, 'b-', 'LineWidth', 1.8, 'DisplayName', 'EMSD (optimized)');
semilogx(Om_vec(ok_base), TF_base_dB, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Baseline (no circuit)');
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('T_F (dB)', 'FontName', fontName, 'FontSize', fsLab);
legend('Location', 'best', 'FontName', fontName, 'FontSize', 10);
title('Force Transmissibility (dB)', 'FontName', fontName, 'FontSize', fsTit);
yline(0, 'k--', 'LineWidth', 1.0);

% 标注峰值抑制
[peak_opt, idx_peak] = max(TF_opt(ok_opt));
[peak_base, ~] = max(TF_base(ok_base));
peak_opt_dB = 20*log10(peak_opt);
plot(Om_vec(find(ok_opt,1)+idx_peak-1), peak_opt_dB, 'bo', 'MarkerSize', 8, 'LineWidth', 1.5);
text(Om_vec(find(ok_opt,1)+idx_peak-1)*1.3, peak_opt_dB+1, ...
    sprintf('Peak=%.1f dB', peak_opt_dB), 'FontName', fontName, 'FontSize', 9, 'Color', 'b');

subplot(1,2,2);
hold on; box on; grid on;
plot(Om_vec(ok_opt), maxMu_opt(ok_opt), 'b.', 'MarkerSize', 6);
plot(Om_vec(ok_base), maxMu_base(ok_base), 'r.', 'MarkerSize', 6);
yline(1.002, 'k--', 'LineWidth', 1.5);
yline(1.0, ':', 'Color', [0.4 0.4 0.4]);
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('max|\mu|', 'FontName', fontName, 'FontSize', fsLab);
legend({'EMSD', 'Baseline', 'tol=1.002', '|\mu|=1'}, ...
    'Location', 'best', 'FontName', fontName, 'FontSize', 9);
title('Floquet Stability', 'FontName', fontName, 'FontSize', fsTit);
set(gca, 'XScale', 'log');
ylim([0.6, 1.8]);

sgtitle(sprintf(['Force Transmissibility & Stability: EMSD vs Baseline\n' ...
    'TF_{peak}: %.1f dB \\rightarrow %.1f dB (%.1f%% reduction),  F_w=%.4f'], ...
    20*log10(peak_base), 20*log10(peak_opt), (1-peak_opt/peak_base)*100, Fw_fixed), ...
    'FontName', fontName, 'FontSize', 14);

if export_pdf
    exportgraphics(gcf, fullfile(out_dir, 'Fig5_TF_comparison.pdf'), 'ContentType', 'vector');
    fprintf('  -> Fig5_TF_comparison.pdf (peak reduction: %.1f%%)\n', (1-peak_opt/peak_base)*100);
end

%% ================================================================
%%  图3: 电路参数扫描 FRF + Floquet —— 第5.2节
%% ================================================================
fprintf('--- 图3: 电路参数扫描 ---\n');

Nvals = 18;
Nw_scan = 120;
Om_scan = logspace(log10(Om_min), log10(Om_max), Nw_scan).';
mult_vec = logspace(log10(0.2), log10(5.0), Nvals);
Nt_floquet = 600;
tol_stable = 1.002;

param_configs = {
    mult_vec * sigma_opt,  '\sigma',   10;
    mult_vec * kap_e_opt,  '\kappa_e', 8;
    mult_vec * kap_c_opt,  '\kappa_c', 9;
};

results = cell(3,1);

for ip = 1:3
    param_vals = param_configs{ip, 1};
    param_name = param_configs{ip, 2};
    param_idx  = param_configs{ip, 3};

    TF_peak  = nan(Nvals, 1);
    maxMu_pk = nan(Nvals, 1);
    stable_pct = nan(Nvals, 1);
    TF_all   = cell(Nvals, 1);
    maxMu_all = cell(Nvals, 1);

    fprintf('  扫描 %s (%d 个值)...\n', param_name, Nvals);

    for iv = 1:Nvals
        sysP = sysP_opt;
        sysP(param_idx) = param_vals(iv);

        TF = nan(Nw_scan, 1);
        maxMu = nan(Nw_scan, 1);
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
            maxMu(j) = compute_floquet_fast(xc, sysP, Om, Nt_floquet);
        end

        TF_all{iv} = TF;
        maxMu_all{iv} = maxMu;
        if nnz(ok) >= 5
            TF_peak(iv)  = max(TF(ok));
            maxMu_pk(iv) = max(maxMu(ok));
            stable_pct(iv) = 100 * sum(maxMu(ok) < tol_stable) / nnz(ok);
        end
    end

    results{ip} = struct('param_vals', param_vals, 'param_name', param_name, ...
        'TF_peak', TF_peak, 'maxMu_pk', maxMu_pk, 'stable_pct', stable_pct, ...
        'TF_all', {TF_all}, 'maxMu_all', {maxMu_all}, 'Om_scan', Om_scan);
end

% 出图
figure('Color','w','Position',[40 40 1350 780]);
for ip = 1:3
    r = results{ip};
    vals = r.param_vals;
    name = r.param_name;
    ref_val = [sigma_opt, kap_e_opt, kap_c_opt];
    ref_val = ref_val(ip);

    show_idx = round(linspace(1, Nvals, 6));
    cmap = lines(6);

    subplot(3, 2, 2*ip-1);
    hold on; box on; grid on;
    for k = 1:6
        idx = show_idx(k);
        tf = r.TF_all{idx};
        ok_pts = isfinite(tf);
        if any(ok_pts)
            semilogx(Om_scan(ok_pts), 20*log10(max(tf(ok_pts), 1e-12)), 'Color', cmap(k,:), 'LineWidth', 1.3, ...
                'DisplayName', sprintf('%s=%.3f', name, vals(idx)));
        end
    end
    xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
    ylabel('T_F (dB)', 'FontName', fontName, 'FontSize', fsLab);
    title(sprintf('FRF vs %s (ref=%.4f)', name, ref_val), ...
        'FontName', fontName, 'FontSize', fsTit);
    legend('Location', 'best', 'FontName', fontName, 'FontSize', 8);

    subplot(3, 2, 2*ip);
    hold on; box on; grid on;
    for k = 1:6
        idx = show_idx(k);
        mu = r.maxMu_all{idx};
        ok_pts = isfinite(mu);
        if any(ok_pts)
            semilogx(Om_scan(ok_pts), mu(ok_pts), 'Color', cmap(k,:), 'LineWidth', 1.2, ...
                'DisplayName', sprintf('%s=%.3f', name, vals(idx)));
        end
    end
    yline(tol_stable, 'k--', 'LineWidth', 1.2);
    yline(1.0, ':', 'Color', [0.4 0.4 0.4]);
    xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
    ylabel('max|\mu|', 'FontName', fontName, 'FontSize', fsLab);
    title(sprintf('Floquet vs %s', name), 'FontName', fontName, 'FontSize', fsTit);
    ylim([0.5, 5.5]);
end
sgtitle(sprintf(['Circuit Parameter Sensitivity\n' ...
    'Ref: \\sigma=%.4f, \\kappa_e=%.4f, \\kappa_c=%.4f'], sigma_opt, kap_e_opt, kap_c_opt), ...
    'FontName', fontName, 'FontSize', 14);

if export_pdf
    exportgraphics(gcf, fullfile(out_dir, 'Fig5_Circuit_Param_FRF.pdf'), 'ContentType', 'vector');
    fprintf('  -> Fig5_Circuit_Param_FRF.pdf\n');
end

%% 灵敏度汇总图
figure('Color','w','Position',[60 60 1250 420]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
p_names = {'\sigma', '\kappa_e', '\kappa_c'};
refs = [sigma_opt, kap_e_opt, kap_c_opt];

for ip = 1:3
    r = results{ip};
    vals = r.param_vals;
    nexttile;
    yyaxis left;
    semilogx(vals, r.TF_peak, 'b-o', 'LineWidth', 1.3, 'MarkerSize', 4);
    xline(refs(ip), 'b--', 'LineWidth', 1.0);
    ylabel('TF_{peak}', 'FontName', fontName, 'FontSize', fsLab);
    grid on; box on;

    yyaxis right;
    semilogx(vals, r.stable_pct, 'r-s', 'LineWidth', 1.3, 'MarkerSize', 4);
    ylabel('Stable %', 'FontName', fontName, 'FontSize', fsLab);
    ylim([0, 105]);

    xlabel(p_names{ip}, 'FontName', fontName, 'FontSize', fsLab);
    title(sprintf('Sensitivity to %s', p_names{ip}), ...
        'FontName', fontName, 'FontSize', fsTit);
    legend({'TF_{peak}', 'Stable %'}, 'Location', 'best', ...
        'FontName', fontName, 'FontSize', 9);
end
sgtitle('Circuit Parameter Sensitivity Summary', 'FontName', fontName, 'FontSize', 14);

if export_pdf
    exportgraphics(gcf, fullfile(out_dir, 'Fig5_Circuit_Param_Sensitivity.pdf'), 'ContentType', 'vector');
    fprintf('  -> Fig5_Circuit_Param_Sensitivity.pdf\n');
end

%% ================================================================
%%  图4: (Fw, Omega) 稳定性边界图 —— 第6.2节
%% ================================================================
fprintf('--- 图4: 稳定性边界图 ---\n');

N_Fw   = 20;
N_Om   = 30;
Fw_vec = logspace(log10(0.001), log10(0.05), N_Fw);
Om_map = logspace(log10(0.2), log10(6.0), N_Om);

fprintf('  计算优化参数稳定性...\n');
[TF_map1, maxMu_map1, ok_map1] = compute_stability_map(...
    sysP_opt, Fw_vec, Om_map, 600);

fprintf('  计算基线稳定性...\n');
[TF_map2, maxMu_map2, ok_map2] = compute_stability_map(...
    sysP_base, Fw_vec, Om_map, 600);

[FF, OO] = meshgrid(Fw_vec, Om_map);
FF = FF'; OO = OO';

figure('Color','w','Position',[50 50 1400 560]);

subplot(2,3,1);
pcolor(FF, OO, log10(max(min(TF_map1, 100), 0.01)));
shading flat; colorbar;
set(gca, 'XScale', 'log', 'YScale', 'log');
colormap(gca, jet);
xlabel('F_w', 'FontName', fontName, 'FontSize', fsLab);
ylabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
title('log_{10} TF (EMSD)', 'FontName', fontName, 'FontSize', fsTit);

subplot(2,3,2);
maxMu_plot1 = maxMu_map1; maxMu_plot1(~ok_map1) = NaN;
pcolor(FF, OO, maxMu_plot1);
shading flat; colorbar; caxis([0.8, 1.3]);
set(gca, 'XScale', 'log', 'YScale', 'log');
hold on;
contour(FF, OO, maxMu_map1 <= tol_stable & ok_map1, [0.5, 0.5], 'k-', 'LineWidth', 1.5);
contour(FF, OO, maxMu_map1 <= 1.0 & ok_map1, [0.5, 0.5], 'k--', 'LineWidth', 1.0);
hold off;
xlabel('F_w', 'FontName', fontName, 'FontSize', fsLab);
ylabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
title('Floquet max|\mu| (EMSD)', 'FontName', fontName, 'FontSize', fsTit);

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

subplot(2,3,4);
pcolor(FF, OO, log10(max(min(TF_map2, 100), 0.01)));
shading flat; colorbar;
set(gca, 'XScale', 'log', 'YScale', 'log');
colormap(gca, jet);
xlabel('F_w', 'FontName', fontName, 'FontSize', fsLab);
ylabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
title('log_{10} TF (Baseline)', 'FontName', fontName, 'FontSize', fsTit);

subplot(2,3,5);
maxMu_plot2 = maxMu_map2; maxMu_plot2(~ok_map2) = NaN;
pcolor(FF, OO, maxMu_plot2);
shading flat; colorbar; caxis([0.8, 1.3]);
set(gca, 'XScale', 'log', 'YScale', 'log');
hold on;
contour(FF, OO, maxMu_map2 <= tol_stable & ok_map2, [0.5, 0.5], 'k-', 'LineWidth', 1.5);
contour(FF, OO, maxMu_map2 <= 1.0 & ok_map2, [0.5, 0.5], 'k--', 'LineWidth', 1.0);
hold off;
xlabel('F_w', 'FontName', fontName, 'FontSize', fsLab);
ylabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
title('Floquet max|\mu| (Baseline)', 'FontName', fontName, 'FontSize', fsTit);

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

sgtitle(sprintf(['Stability Boundary Map in (F_w, \\Omega) Plane\n' ...
    'EMSD: \\sigma=%.4f, \\kappa_e=%.4f, \\kappa_c=%.4f'], ...
    sigma_opt, kap_e_opt, kap_c_opt), 'FontName', fontName, 'FontSize', 14);

if export_pdf
    exportgraphics(gcf, fullfile(out_dir, 'Fig6_Stability_Boundary_Map.pdf'), 'ContentType', 'vector');
    fprintf('  -> Fig6_Stability_Boundary_Map.pdf\n');
end

%% ================================================================
%%  图5: 定频扫力 + Floquet —— 第6.1节
%% ================================================================
fprintf('--- 图5: 定频扫力 + Floquet ---\n');

% 选3个代表性频率
Om_demo = [0.5, 1.0, 2.0];
N_Fw_demo = 30;
Fw_demo = logspace(log10(0.001), log10(0.05), N_Fw_demo);
Nt_fq = 600;

figure('Color','w','Position',[50 50 1400 420]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

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
        % 响应幅值（x1的一次谐波）
        A_resp(iFw) = hypot(xc(2), xc(3));
        maxMu_arr(iFw) = compute_floquet_fast(xc, sysP_opt, Om, Nt_fq);
        stable_arr(iFw) = maxMu_arr(iFw) < tol_stable;
    end

    nexttile;
    yyaxis left;
    ok_pts = isfinite(A_resp) & A_resp > 0;
    scatter(Fw_demo(stable_arr & ok_pts), A_resp(stable_arr & ok_pts), ...
        20, 'b', 'filled', 'DisplayName', 'Stable');
    hold on;
    scatter(Fw_demo(~stable_arr & ok_pts), A_resp(~stable_arr & ok_pts), ...
        20, 'r', 'DisplayName', 'Unstable');
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('F_w', 'FontName', fontName, 'FontSize', fsLab);
    ylabel('|x_1|_{amp}', 'FontName', fontName, 'FontSize', fsLab);
    grid on; box on;
    title(sprintf('\\Omega=%.2f', Om), 'FontName', fontName, 'FontSize', fsTit);

    yyaxis right;
    plot(Fw_demo(ok_pts), maxMu_arr(ok_pts), 'k.-', 'MarkerSize', 6, 'LineWidth', 0.8);
    yline(1.002, 'k--', 'LineWidth', 1.0);
    ylabel('max|\mu|', 'FontName', fontName, 'FontSize', fsLab);
    legend('Location', 'best', 'FontName', fontName, 'FontSize', 7);
    ylim([0.5, 2.5]);
end

sgtitle('Force Sweep: Response Amplitude & Floquet Stability', ...
    'FontName', fontName, 'FontSize', 14);

if export_pdf
    exportgraphics(gcf, fullfile(out_dir, 'Fig6_Force_Sweep_Floquet.pdf'), 'ContentType', 'vector');
    fprintf('  -> Fig6_Force_Sweep_Floquet.pdf\n');
end

%% ================================================================
%%  完成
%% ================================================================
fprintf('\n========================================\n');
fprintf('  全部图片已保存到: %s\n', out_dir);
fprintf('========================================\n');

% 列出所有生成的文件
files = dir(fullfile(out_dir, '*.pdf'));
fprintf('\n生成的图片文件 (%d 个):\n', length(files));
for k = 1:length(files)
    fprintf('  [%d] %s  (%.1f KB)\n', k, files(k).name, files(k).bytes/1024);
end

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
