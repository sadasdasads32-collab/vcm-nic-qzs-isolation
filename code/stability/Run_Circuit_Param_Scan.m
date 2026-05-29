%% Run_Circuit_Param_Scan.m
% =========================================================================
% 电路参数对稳定性的影响：扫描 sigma, kap_e, kap_c
%
% 以优化参数为中心，每个参数单独扫描 18 个值（保持其他两个固定），
% 逐频点 HBM → FRF + Floquet，分析 TF 峰值和稳定性的参数敏感性
% =========================================================================

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'lib'));

%% ===================== 参考点（优化参数）=====================

sigma_ref = 1.1506;
kap_e_ref = 1.5222;
kap_c_ref = 0.5743;

%% ===================== 机械参数 =====================

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

sysP0 = [be1, be2, mu, al1, gamma1, ze1, lam_phys, ...
         kap_e_ref, kap_c_ref, sigma_ref, gamma2];

%% ===================== 扫描设置 =====================

Fw_fixed = 0.008;
Om_min   = 0.2;
Om_max   = 6.0;
Nw_scan  = 120;  % 频点（比优化时略稀疏以加速）

Nvals    = 18;
Nt_floquet = 600;
tol_stable = 1.002;

global FixedOmega Fw
Fw = Fw_fixed;
FixedOmega = [];

Om_scan = logspace(log10(Om_min), log10(Om_max), Nw_scan).';

% 扫描倍数范围：0.2x ~ 5x 参考值
mult_vec = logspace(log10(0.2), log10(5.0), Nvals);

% 三个参数的扫描值和标签
param_configs = {
    mult_vec * sigma_ref,  '\sigma',   10;   % sigma扫描
    mult_vec * kap_e_ref,  '\kappa_e', 8;    % kap_e扫描
    mult_vec * kap_c_ref,  '\kappa_c', 9;    % kap_c扫描
};

results = cell(3, 1);

fprintf('========== 电路参数敏感性扫描 ==========\n');
fprintf('参考值: sigma=%.4f, kap_e=%.4f, kap_c=%.4f\n', ...
    sigma_ref, kap_e_ref, kap_c_ref);
fprintf('每参数 %d 个值 x %d 频点\n\n', Nvals, Nw_scan);

%% ===================== 逐参数扫描 =====================

for ip = 1:3
    param_vals = param_configs{ip, 1};
    param_name = param_configs{ip, 2};
    param_idx  = param_configs{ip, 3};  % 在 sysP 中的位置

    TF_peak  = nan(Nvals, 1);
    maxMu_pk = nan(Nvals, 1);
    stable_pct = nan(Nvals, 1);
    TF_all   = cell(Nvals, 1);
    maxMu_all = cell(Nvals, 1);

    fprintf('--- 扫描 %s (%d 个值) ---\n', param_name, Nvals);

    for iv = 1:Nvals
        sysP = sysP0;
        sysP(param_idx) = param_vals(iv);

        % 快速 FRF 扫频
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

        if mod(iv, 6) == 0 || iv == 1
            fprintf('  %s=%7.4f (%2d/%d): TFpk=%.3f, max|mu|=%.3f, stb=%.0f%%\n', ...
                param_name, param_vals(iv), iv, Nvals, TF_peak(iv), maxMu_pk(iv), stable_pct(iv));
        end
    end

    results{ip} = struct('param_vals', param_vals, 'param_name', param_name, ...
        'TF_peak', TF_peak, 'maxMu_pk', maxMu_pk, 'stable_pct', stable_pct, ...
        'TF_all', {TF_all}, 'maxMu_all', {maxMu_all}, 'Om_scan', Om_scan);
end

%% ===================== 出图 =====================

fontName = 'Times New Roman';
fsLab = 11; fsTit = 12;

figure('Color','w','Position',[40 40 1350 780]);

for ip = 1:3
    r = results{ip};
    vals = r.param_vals;
    name = r.param_name;
    ref_val = [sigma_ref, kap_e_ref, kap_c_ref];
    ref_val = ref_val(ip);

    % 选 6 条代表性曲线（均匀分布）
    show_idx = round(linspace(1, Nvals, 6));
    cmap = lines(6);

    % --- 左列: TF FRF ---
    subplot(3, 2, 2*ip-1);
    hold on; box on; grid on;
    for k = 1:6
        idx = show_idx(k);
        tf = r.TF_all{idx};
        ok_pts = isfinite(tf);
        if any(ok_pts)
            loglog(Om_scan(ok_pts), tf(ok_pts), 'Color', cmap(k,:), ...
                'LineWidth', 1.3, ...
                'DisplayName', sprintf('%s=%.3f', name, vals(idx)));
        end
    end
    xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
    ylabel('T_F(\Omega)', 'FontName', fontName, 'FontSize', fsLab);
    title(sprintf('FRF vs %s (ref=%.4f)', name, ref_val), ...
        'FontName', fontName, 'FontSize', fsTit);
    legend('Location', 'best', 'FontName', fontName, 'FontSize', 8);

    % --- 右列: Floquet ---
    subplot(3, 2, 2*ip);
    hold on; box on; grid on;
    for k = 1:6
        idx = show_idx(k);
        mu = r.maxMu_all{idx};
        ok_pts = isfinite(mu);
        if any(ok_pts)
            semilogx(Om_scan(ok_pts), mu(ok_pts), 'Color', cmap(k,:), ...
                'LineWidth', 1.2, ...
                'DisplayName', sprintf('%s=%.3f', name, vals(idx)));
        end
    end
    yline(tol_stable, 'k--', 'LineWidth', 1.2);
    yline(1.0, ':', 'Color', [0.4 0.4 0.4]);
    xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
    ylabel('max|\mu|', 'FontName', fontName, 'FontSize', fsLab);
    title(sprintf('Floquet vs %s', name), 'FontName', fontName, 'FontSize', fsTit);
    ylim([0.5, 2.5]);
end

sgtitle(sprintf(['Circuit Parameter Sensitivity\n' ...
    'Ref: \\sigma=%.4f, \\kappa_e=%.4f, \\kappa_c=%.4f, F_w=%.4f'], ...
    sigma_ref, kap_e_ref, kap_c_ref, Fw_fixed), ...
    'FontName', fontName, 'FontSize', 14);

%% ===================== 汇总图: TF峰值 + 稳定性 vs 参数 =====================

figure('Color','w','Position',[60 60 1250 420]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

p_names = {'\sigma', '\kappa_e', '\kappa_c'};
refs = [sigma_ref, kap_e_ref, kap_c_ref];

for ip = 1:3
    r = results{ip};
    vals = r.param_vals;
    nn = norm(vals - min(vals)) / norm(max(vals) - min(vals));

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

sgtitle('Circuit Parameter Sensitivity Summary', ...
    'FontName', fontName, 'FontSize', 14);

%% ===================== 保存 =====================

out_dir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

save(fullfile(out_dir, ['circuit_param_scan_' timestamp '.mat']), ...
    'results', 'sigma_ref', 'kap_e_ref', 'kap_c_ref', ...
    'Fw_fixed', 'Om_scan', 'Nt_floquet', 'tol_stable');

figs = findall(0, 'Type', 'figure');
for k = 1:length(figs)
    fname = sprintf('param_scan_%02d_%s', k, timestamp);
    exportgraphics(figs(k), fullfile(out_dir, [fname '.pdf']), ...
        'ContentType', 'vector');
    fprintf('  导出: %s.pdf\n', fname);
end

fprintf('\n数据已保存到: %s\n', out_dir);
fprintf('图片已导出\n');
