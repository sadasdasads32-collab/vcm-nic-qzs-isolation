%% Generate_Chapter4_Figures.m
% =========================================================================
% 生成第4章"复动力学算子的频率塑形机理分析"所需全部图片
% 使用优化参数: sigma=1.1506, kap_e=1.5222, kap_c=0.5743
% 保存到 e:\项目1\论文图\
% =========================================================================

clc; clear; close all;
init_path();

out_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'output', 'journal_figures');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

fontName = 'Times New Roman';
fsLab = 12; fsTit = 13;

%% ==================== 优化参数 ====================

sigma_opt = 1.1506;
kap_e_opt = 1.5222;
kap_c_opt = 0.5743;
theta     = sqrt(0.18);  % lam_phys = 0.18

fprintf('========== Chapter 4: Complex Operator Figures ==========\n');
fprintf('Parameters: sigma=%.4f, kap_e=%.4f, kap_c=%.4f, theta=%.4f\n\n', ...
    sigma_opt, kap_e_opt, kap_c_opt, theta);

%% ==================== 频率范围 ====================

Om_min = 0.05;
Om_max = 10.0;
N_om = 2000;
Omega = logspace(log10(Om_min), log10(Om_max), N_om).';

%% ==================== 计算 K(Omega) ====================

Delta = kap_e_opt * Omega.^2 - kap_c_opt;
Den   = Delta.^2 + (sigma_opt * Omega).^2;

Kr = theta^2 .* Omega.^2 .* Delta ./ Den;   % K 实部
Ki = theta^2 .* sigma_opt .* Omega.^3 ./ Den;  % K 虚部

Omega_e = sqrt(kap_c_opt / kap_e_opt);  % 电路特征频率

% 等效参数
C_eq = Ki ./ Omega;           % Ceq = Im[K]/Omega (damping)
M_eq = -Kr ./ (Omega.^2);     % Meq = -Re[K]/Omega^2 (inertia, negative sign for virtual inertia)
K_eq_stiff = Kr;              % Keq = Re[K] (stiffness contribution)

fprintf('Circuit characteristic frequency: Omega_e = %.4f\n', Omega_e);

% 渐近分析
M_eq_low = theta^2 / kap_c_opt;   % 低频虚拟惯性
Kr_high  = theta^2 / kap_e_opt;   % 高频常数极限

fprintf('Low-freq virtual inertia: M_eq(0) = %.4f\n', M_eq_low);
fprintf('High-freq limit: K(inf) -> %.4f\n', Kr_high);

%% ================================================================
%%  图4.1: K_R(Omega) 和 K_I(Omega) 及频带分区
%% ================================================================
fprintf('--- Fig4_1: K_R & K_I with frequency regions ---\n');

figure('Color','w','Position',[100 100 1000 500]);
hold on; grid on; box on;
set(gca, 'XScale', 'log');

plot(Omega, Kr, 'b-', 'LineWidth', 2.2);
plot(Omega, Ki, 'r-', 'LineWidth', 2.2);
xline(Omega_e, 'k--', 'LineWidth', 1.5);
yline(0, 'k-', 'LineWidth', 0.8);

% 标注三个频段区域
yl = ylim();

% 虚拟惯性区（低频）
xf1 = Om_min * 1.5;
xf2 = Omega_e * 0.55;
fill_x = [xf1 xf2 xf2 xf1];
fill_y = [yl(1) yl(1) yl(2) yl(2)];
patch(fill_x, fill_y, [0.6 0.8 1.0], 'FaceAlpha', 0.12, 'EdgeColor', 'none');
text(sqrt(xf1*xf2), yl(2)*0.88, 'Virtual Inertia', ...
    'FontName', fontName, 'FontSize', 10, 'Color', 'b', ...
    'HorizontalAlignment', 'center');

% 阻尼塑形区（中频）
xf3 = Omega_e * 0.65;
xf4 = Omega_e * 2.8;
fill_x = [xf3 xf4 xf4 xf3];
fill_y = [yl(1) yl(1) yl(2) yl(2)];
patch(fill_x, fill_y, [1.0 0.7 0.7], 'FaceAlpha', 0.12, 'EdgeColor', 'none');
text(sqrt(xf3*xf4), yl(2)*0.82, 'Damping Shaping', ...
    'FontName', fontName, 'FontSize', 10, 'Color', 'r', ...
    'HorizontalAlignment', 'center');

% 类刚度区（高频）
xf5 = Omega_e * 2.3;
xf6 = Omega_e * 7;
fill_x = [xf5 xf6 xf6 xf5];
fill_y = [yl(1) yl(1) yl(2) yl(2)];
patch(fill_x, fill_y, [0.7 1.0 0.7], 'FaceAlpha', 0.12, 'EdgeColor', 'none');
text(sqrt(xf5*xf6), yl(2)*0.88, 'Stiffness-like', ...
    'FontName', fontName, 'FontSize', 10, 'Color', [0 0.5 0], ...
    'HorizontalAlignment', 'center');

xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('K(\Omega)', 'FontName', fontName, 'FontSize', fsLab);
title(sprintf(['Real and Imaginary Parts of Complex Operator K(\\Omega)\n' ...
    '\\sigma=%.4f, \\kappa_e=%.4f, \\kappa_c=%.4f, \\Omega_e=%.3f'], ...
    sigma_opt, kap_e_opt, kap_c_opt, Omega_e), ...
    'FontName', fontName, 'FontSize', fsTit);

legend({'K_R(\Omega)', 'K_I(\Omega)', '\Omega_e'}, ...
    'Location', 'northeast', 'FontName', fontName, 'FontSize', 11);

xlim([Om_min Om_max]);

exportgraphics(gcf, fullfile(out_dir, 'Fig4_1_Kr_Ki.pdf'), 'ContentType', 'vector');
fprintf('  -> Fig4_1_Kr_Ki.pdf\n');

%% ================================================================
%%  图4.2: 等效惯性 Meq、等效阻尼 Ceq、等效刚度 Keq 三子图
%% ================================================================
fprintf('--- Fig4_2: Equivalent parameters Meq, Ceq, Keq ---\n');

figure('Color','w','Position',[100 100 1000 800]);

% (a) 等效惯性
subplot(3,1,1);
hold on; grid on; box on;
set(gca, 'XScale', 'log');

% 正虚拟惯性区（低频：Kr<0, Meq>0）
pos_iner = M_eq > 0;
plot(Omega(pos_iner), M_eq(pos_iner), 'b-', 'LineWidth', 1.8);
% 负区（高频：Kr>0, Meq<0, 对应正刚度贡献）
neg_iner = M_eq <= 0;
if any(neg_iner)
    plot(Omega(neg_iner), M_eq(neg_iner), 'b--', 'LineWidth', 1.2);
end
xline(Omega_e, 'k--', 'LineWidth', 1.5);
yline(0, 'k-', 'LineWidth', 0.8);

% 标注
text(Omega_e*0.3, max(M_eq)*0.85, ...
    sprintf('M_{eq}(0) = %.3f', M_eq_low), ...
    'FontName', fontName, 'FontSize', 10, 'Color', 'b');

xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('M_{eq}(\Omega)', 'FontName', fontName, 'FontSize', fsLab);
title('Equivalent Virtual Inertia: M_{eq} = -K_R/\Omega^2', ...
    'FontName', fontName, 'FontSize', fsTit);
legend({'M_{eq}>0 (virtual inertia)', 'M_{eq}<0 (mass loading)'}, ...
    'Location', 'best', 'FontName', fontName, 'FontSize', 9);

% (b) 等效阻尼
subplot(3,1,2);
hold on; grid on; box on;
set(gca, 'XScale', 'log');
plot(Omega, C_eq, 'r-', 'LineWidth', 1.8);
xline(Omega_e, 'k--', 'LineWidth', 1.5);

[Ce_peak, idx_ce] = max(C_eq);
plot(Omega(idx_ce), Ce_peak, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
text(Omega(idx_ce)*1.3, Ce_peak*0.85, ...
    sprintf('C_{eq,max}=%.4f at \\Omega=%.3f', Ce_peak, Omega(idx_ce)), ...
    'FontName', fontName, 'FontSize', 10, 'Color', 'r');

xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('C_{eq}(\Omega)', 'FontName', fontName, 'FontSize', fsLab);
title('Equivalent Damping: C_{eq} = K_I/\Omega', ...
    'FontName', fontName, 'FontSize', fsTit);

% (c) 等效刚度
subplot(3,1,3);
hold on; grid on; box on;
set(gca, 'XScale', 'log');
plot(Omega, K_eq_stiff, 'Color', [0 0.5 0], 'LineWidth', 1.8);
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

xlim([Om_min Om_max]);

exportgraphics(gcf, fullfile(out_dir, 'Fig4_2_Equivalent_Params.pdf'), 'ContentType', 'vector');
fprintf('  -> Fig4_2_Equivalent_Params.pdf\n');

%% ================================================================
%%  图4.3: 归一化 K_R 和 K_I
%% ================================================================
fprintf('--- Fig4_3: Normalized K_R & K_I ---\n');

Kr_norm = Kr / max(abs(Kr));
Ki_norm = Ki / max(abs(Ki));

figure('Color','w','Position',[100 100 900 460]);
hold on; grid on; box on;
set(gca, 'XScale', 'log');

plot(Omega, Kr_norm, 'b-', 'LineWidth', 2.2);
plot(Omega, Ki_norm, 'r-', 'LineWidth', 2.2);
xline(Omega_e, 'k--', 'LineWidth', 1.5);
yline(0, 'k-', 'LineWidth', 0.8);

% 标注区域
yl = ylim();
text(Omega_e*0.3, 0.85, 'Virtual\nInertia', ...
    'FontName', fontName, 'FontSize', 10, 'Color', 'b', ...
    'HorizontalAlignment', 'center');
text(Omega_e*1.3, 0.85, 'Damping\nShaping', ...
    'FontName', fontName, 'FontSize', 10, 'Color', 'r', ...
    'HorizontalAlignment', 'center');
text(Omega_e*4.5, 0.45, 'Stiffness-like\n(constant)', ...
    'FontName', fontName, 'FontSize', 10, 'Color', [0 0.5 0], ...
    'HorizontalAlignment', 'center');

xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('Normalized Value', 'FontName', fontName, 'FontSize', fsLab);
title(sprintf('Normalized K_R(\\Omega) and K_I(\\Omega), \\Omega_e=%.3f', Omega_e), ...
    'FontName', fontName, 'FontSize', fsTit);
legend({'K_R / max|K_R|', 'K_I / max|K_I|', '\Omega_e'}, ...
    'Location', 'best', 'FontName', fontName, 'FontSize', 10);

xlim([Om_min Om_max]);

exportgraphics(gcf, fullfile(out_dir, 'Fig4_3_Normalized_K.pdf'), 'ContentType', 'vector');
fprintf('  -> Fig4_3_Normalized_K.pdf\n');

%% ================================================================
%%  图4.4: 电路参数扫描对 K(Omega) 形状的影响 —— 第4.3节
%% ================================================================
fprintf('--- Fig4_4: Circuit parameter effects on K(Omega) ---\n');

% 扫描参数
sig_vals = [0.3, 0.6, sigma_opt, 2.0, 4.0];
kpe_vals = [0.3, 0.8, kap_e_opt, 2.5, 5.0];
kpc_vals = [0.1, 0.3, kap_c_opt, 1.0, 2.0];

colors = lines(5);

figure('Color','w','Position',[50 50 1400 900]);

% --- (a, b): sigma 扫描 ---
subplot(3,2,1);
hold on; grid on; box on; set(gca, 'XScale', 'log');
leg_str = {};
for k = 1:5
    s = sig_vals(k);
    D_s = kap_e_opt * Omega.^2 - kap_c_opt;
    Den_s = D_s.^2 + (s * Omega).^2;
    Kr_s = theta^2 .* Omega.^2 .* D_s ./ Den_s;
    plot(Omega, Kr_s, '-', 'Color', colors(k,:), 'LineWidth', 1.5);
    leg_str{end+1} = sprintf('\\sigma=%.2f', s);
end
xline(Omega_e, 'k--', 'LineWidth', 1.0);
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('K_R(\Omega)', 'FontName', fontName, 'FontSize', fsLab);
title('Real Part K_R vs \sigma', 'FontName', fontName, 'FontSize', fsTit);
legend(leg_str, 'Location', 'best', 'FontName', fontName, 'FontSize', 8);

subplot(3,2,2);
hold on; grid on; box on; set(gca, 'XScale', 'log');
for k = 1:5
    s = sig_vals(k);
    D_s = kap_e_opt * Omega.^2 - kap_c_opt;
    Den_s = D_s.^2 + (s * Omega).^2;
    Ki_s = theta^2 .* s .* Omega.^3 ./ Den_s;
    plot(Omega, Ki_s, '-', 'Color', colors(k,:), 'LineWidth', 1.5);
end
xline(Omega_e, 'k--', 'LineWidth', 1.0);
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('K_I(\Omega)', 'FontName', fontName, 'FontSize', fsLab);
title('Imaginary Part K_I vs \sigma', 'FontName', fontName, 'FontSize', fsTit);
legend(leg_str, 'Location', 'best', 'FontName', fontName, 'FontSize', 8);

% --- (c, d): kap_e 扫描 ---
subplot(3,2,3);
hold on; grid on; box on; set(gca, 'XScale', 'log');
leg_str = {};
for k = 1:5
    ke = kpe_vals(k);
    D_k = ke * Omega.^2 - kap_c_opt;
    Den_k = D_k.^2 + (sigma_opt * Omega).^2;
    Kr_k = theta^2 .* Omega.^2 .* D_k ./ Den_k;
    plot(Omega, Kr_k, '-', 'Color', colors(k,:), 'LineWidth', 1.5);
    leg_str{end+1} = sprintf('\\kappa_e=%.2f', ke);
end
xline(Omega_e, 'k--', 'LineWidth', 1.0);
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('K_R(\Omega)', 'FontName', fontName, 'FontSize', fsLab);
title('Real Part K_R vs \kappa_e', 'FontName', fontName, 'FontSize', fsTit);
legend(leg_str, 'Location', 'best', 'FontName', fontName, 'FontSize', 8);

subplot(3,2,4);
hold on; grid on; box on; set(gca, 'XScale', 'log');
for k = 1:5
    ke = kpe_vals(k);
    D_k = ke * Omega.^2 - kap_c_opt;
    Den_k = D_k.^2 + (sigma_opt * Omega).^2;
    Ki_k = theta^2 .* sigma_opt .* Omega.^3 ./ Den_k;
    plot(Omega, Ki_k, '-', 'Color', colors(k,:), 'LineWidth', 1.5);
end
xline(Omega_e, 'k--', 'LineWidth', 1.0);
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('K_I(\Omega)', 'FontName', fontName, 'FontSize', fsLab);
title('Imaginary Part K_I vs \kappa_e', 'FontName', fontName, 'FontSize', fsTit);
legend(leg_str, 'Location', 'best', 'FontName', fontName, 'FontSize', 8);

% --- (e, f): kap_c 扫描 ---
subplot(3,2,5);
hold on; grid on; box on; set(gca, 'XScale', 'log');
leg_str = {};
for k = 1:5
    kc = kpc_vals(k);
    D_c = kap_e_opt * Omega.^2 - kc;
    Den_c = D_c.^2 + (sigma_opt * Omega).^2;
    Kr_c = theta^2 .* Omega.^2 .* D_c ./ Den_c;
    plot(Omega, Kr_c, '-', 'Color', colors(k,:), 'LineWidth', 1.5);
    leg_str{end+1} = sprintf('\\kappa_c=%.2f', kc);
end
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('K_R(\Omega)', 'FontName', fontName, 'FontSize', fsLab);
title('Real Part K_R vs \kappa_c', 'FontName', fontName, 'FontSize', fsTit);
legend(leg_str, 'Location', 'best', 'FontName', fontName, 'FontSize', 8);

subplot(3,2,6);
hold on; grid on; box on; set(gca, 'XScale', 'log');
for k = 1:5
    kc = kpc_vals(k);
    D_c = kap_e_opt * Omega.^2 - kc;
    Den_c = D_c.^2 + (sigma_opt * Omega).^2;
    Ki_c = theta^2 .* sigma_opt .* Omega.^3 ./ Den_c;
    plot(Omega, Ki_c, '-', 'Color', colors(k,:), 'LineWidth', 1.5);
end
xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('K_I(\Omega)', 'FontName', fontName, 'FontSize', fsLab);
title('Imaginary Part K_I vs \kappa_c', 'FontName', fontName, 'FontSize', fsTit);
legend(leg_str, 'Location', 'best', 'FontName', fontName, 'FontSize', 8);

sgtitle(sprintf(['Effect of Circuit Parameters on K(\\Omega) Shape\n' ...
    'Reference: \\sigma=%.4f, \\kappa_e=%.4f, \\kappa_c=%.4f, \\Omega_e=%.3f'], ...
    sigma_opt, kap_e_opt, kap_c_opt, Omega_e), ...
    'FontName', fontName, 'FontSize', 14);

exportgraphics(gcf, fullfile(out_dir, 'Fig4_4_Parameter_Effects.pdf'), 'ContentType', 'vector');
fprintf('  -> Fig4_4_Parameter_Effects.pdf\n');

%% ================================================================
%%  图4.5: 频带分工与参数映射总结图
%% ================================================================
fprintf('--- Fig4_5: Frequency band mapping summary ---\n');

figure('Color','w','Position',[100 100 1100 480]);

% 左图：等效阻尼峰值位置 vs kap_e, kap_c
subplot(1,2,1);
hold on; grid on; box on;

kap_e_range = linspace(0.5, 4.0, 40);
kap_c_range = linspace(0.1, 2.0, 40);
[KE, KC] = meshgrid(kap_e_range, kap_c_range);
Omega_e_map = sqrt(KC ./ KE);

% 分母零点安全边界: kap_e*Omega^2 - kap_c = 0 => Omega = sqrt(kap_c/kap_e) = Omega_e
% 在 Omega_e 处，如果 sigma 不够大，分母 = (sigma*Omega_e)^2，K 可能很大
% 安全条件：分母在频段 [0.2, 6.0] 内不趋近于零
% 即 kap_e*Omega^2 - kap_c = 0 的解 Omega_0 = sqrt(kap_c/kap_e) 需在 [0.2, 6.0] 之外
% 或 sigma 足够大使得分母数量级不低于某阈值

imagesc(kap_e_range, kap_c_range, Omega_e_map);
colorbar;
colormap(gca, jet);
hold on;
plot(kap_e_opt, kap_c_opt, 'wo', 'MarkerSize', 12, 'LineWidth', 2);
plot(kap_e_opt, kap_c_opt, 'k+', 'MarkerSize', 12, 'LineWidth', 2);
text(kap_e_opt+0.1, kap_c_opt+0.06, 'Optimal', ...
    'FontName', fontName, 'FontSize', 10, 'Color', 'w', 'FontWeight', 'bold');

% 画 Omega_e = 0.2, 0.8, 2.0, 6.0 的等高线
contour(KE, KC, Omega_e_map, [0.2 0.8 2.0 6.0], 'w-', 'LineWidth', 1.2);

xlabel('\kappa_e', 'FontName', fontName, 'FontSize', fsLab);
ylabel('\kappa_c', 'FontName', fontName, 'FontSize', fsLab);
title('\Omega_e = \surd(\kappa_c/\kappa_e) Map', ...
    'FontName', fontName, 'FontSize', fsTit);

% 右图：频带分工示意
subplot(1,2,2);
hold on; grid on; box on;
set(gca, 'XScale', 'log');

% 画归一化 Kr 和 Ki 的简化示意（仅用参考参数）
Kr_demo = Kr / max(abs(Kr));
Ki_demo = Ki / max(abs(Ki));

plot(Omega, Kr_demo, 'b-', 'LineWidth', 2.2);
plot(Omega, Ki_demo, 'r-', 'LineWidth', 2.2);
xline(Omega_e, 'k--', 'LineWidth', 1.5);
yline(0, 'k-', 'LineWidth', 0.8);

% 标注参数控制区域
yl = ylim();

% sigma 控制区
xf_mid = Omega_e * 0.6;
xf_mid2 = Omega_e * 3.0;
patch([xf_mid xf_mid2 xf_mid2 xf_mid], [yl(1) yl(1) yl(2) yl(2)], ...
    [1 0.8 0.8], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
text(Omega_e*0.95, yl(2)*0.9, '\sigma: damping peak', ...
    'FontName', fontName, 'FontSize', 10, 'Color', 'r', 'FontWeight', 'bold');

% kappa_c 控制区 (低频)
patch([Om_min Om_min Omega_e*0.4 Omega_e*0.4], ...
    [yl(1) yl(2) yl(2) yl(1)], ...
    [0.6 0.8 1.0], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
text(Om_min*2, yl(2)*0.75, '\kappa_c:\nvirtual inertia', ...
    'FontName', fontName, 'FontSize', 9, 'Color', 'b');

% kappa_e 控制区 (全频+高频)
patch([Omega_e*0.5 Omega_e*6 Omega_e*6 Omega_e*0.5], ...
    [yl(1) yl(1) 0 0], ...
    [0.7 1.0 0.7], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
text(Omega_e*2.5, -0.15, '\kappa_e: transition & high-freq limit', ...
    'FontName', fontName, 'FontSize', 9, 'Color', [0 0.5 0]);

xlabel('\Omega', 'FontName', fontName, 'FontSize', fsLab);
ylabel('Normalized K(\Omega)', 'FontName', fontName, 'FontSize', fsLab);
title('Frequency-Band Parameter-to-Effect Mapping', ...
    'FontName', fontName, 'FontSize', fsTit);
legend({'K_R (norm)', 'K_I (norm)', '\Omega_e'}, ...
    'Location', 'southwest', 'FontName', fontName, 'FontSize', 9);

sgtitle(sprintf(['Frequency Shaping via K(\\Omega): Band Assignment & Parameter Control\n' ...
    '\\sigma \\rightarrow Damping,  \\kappa_e \\rightarrow Transition/High-freq,  ' ...
    '\\kappa_c \\rightarrow Low-freq Inertia'], ...
    sigma_opt, kap_e_opt, kap_c_opt), ...
    'FontName', fontName, 'FontSize', 13);

exportgraphics(gcf, fullfile(out_dir, 'Fig4_5_FreqBand_Mapping.pdf'), 'ContentType', 'vector');
fprintf('  -> Fig4_5_FreqBand_Mapping.pdf\n');

%% ================================================================
%%  完成
%% ================================================================
fprintf('\n========== All Chapter 4 figures saved ==========\n');

files = dir(fullfile(out_dir, 'Fig4_*.pdf'));
fprintf('\nGenerated Chapter 4 figures (%d files):\n', length(files));
total_size = 0;
for k = 1:length(files)
    fprintf('  [%d] %s  (%.1f KB)\n', k, files(k).name, files(k).bytes/1024);
    total_size = total_size + files(k).bytes;
end
fprintf('  Total: %.1f KB\n', total_size/1024);

% Also list all figures in the directory
all_files = dir(fullfile(out_dir, '*.pdf'));
fprintf('\nAll figures in %s (%d files):\n', out_dir, length(all_files));
for k = 1:length(all_files)
    fprintf('  %s\n', all_files(k).name);
end
