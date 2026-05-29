%% Run_Bifurcation_Classification.m
% =========================================================================
% Bifurcation type classification for the VCM-NIC QZS system
%
% Runs force sweep at 3 representative frequencies (Omega = 0.5, 1.0, 2.0)
% and classifies the bifurcation type at each stability boundary using
% Floquet multiplier analysis.
%
% Bifurcation types:
%   Fold (Saddle-node):  Real multiplier exits at (+1, 0)
%   Flip (Period-doubling): Real multiplier exits at (-1, 0)
%   Neimark-Sacker:       Complex conjugate pair exits unit circle
%
% Outputs:
%   1. Annotated A-F figure with bifurcation type labels at stability boundaries
%   2. Supplementary Argand diagram of dominant multiplier trajectory
% =========================================================================

clc; clear; close all;
init_path();

%% ==================== Parameters ====================
sigma_opt = 1.1506;
kap_e_opt = 1.5222;
kap_c_opt = 0.5743;
lam_phys  = 0.18;

mu   = 0.2;   beta_m = 2.0;   K1 = 1.0;   K2 = 0.2;
U    = 2.0;   Lg = 4/9;       v  = 2.5;
alpha1 = v    - 2*K1*(1-Lg)/Lg;
alpha2 = beta_m - 2*K2*(1-Lg)/Lg;
gamma1_val = K1 / (U^2 * Lg^3);
gamma2_val = K2 / (U^2 * Lg^3);
theta = sqrt(max(lam_phys, 0));
ze1 = 0.05;
be1 = 1.0;
al1 = alpha1 - be1;
be2 = alpha2;

sysP_opt = [be1, be2, mu, al1, gamma1_val, ze1, lam_phys, ...
            kap_e_opt, kap_c_opt, sigma_opt, gamma2_val];

tol_stable = 1.002;
fontName = 'Times New Roman';
fsLab = 12; fsTit = 13;

global FixedOmega Fw

fprintf('========================================\n');
fprintf('  Bifurcation Classification Analysis\n');
fprintf('========================================\n\n');

%% ==================== Force Sweep at 3 Frequencies ====================
Om_demo = [0.5, 1.0, 2.0];
N_Fw = 30;
Fw_vec = logspace(log10(0.001), log10(0.05), N_Fw);
Nt_floquet = 600;

% Store results for figure annotation
bif_data = cell(3, 1);
mu_trajectories = cell(3, 1);

for iOm = 1:3
    Om = Om_demo(iOm);
    FixedOmega = Om;
    Fw = Fw_vec(1);

    fprintf('--- Omega = %.2f ---\n', Om);

    A_resp     = nan(N_Fw, 1);
    maxMu_arr  = nan(N_Fw, 1);
    stable_arr = false(N_Fw, 1);
    bif_types  = cell(N_Fw, 1);
    mu_all_store = cell(N_Fw, 1);

    y_guess = [zeros(15,1); Fw_vec(1)];
    fail_count = 0;
    prev_mu_all = [];  % Store previous step's full multiplier set

    for iFw = 1:N_Fw
        Fw = Fw_vec(iFw);
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

        % Get full Floquet multiplier set (all 6 eigenvalues)
        mu_all = get_full_multipliers_fast(xc, sysP_opt, Om, Nt_floquet);
        mu_all_store{iFw} = mu_all;
        maxMu_arr(iFw) = max(abs(mu_all));
        stable_arr(iFw) = maxMu_arr(iFw) < tol_stable;

        % Classify bifurcation if we have two consecutive steps
        if iFw > 1 && ~isempty(prev_mu_all)
            [bif_type, detail] = classify_bifurcation(prev_mu_all, mu_all, tol_stable);
            bif_types{iFw} = bif_type;
            if ~strcmp(bif_type, 'none')
                fprintf('  Fw=%.4f -> %.4f: %s bifurcation detected (max|mu|=%.4f)\n', ...
                    Fw_vec(iFw-1), Fw_vec(iFw), bif_type, max(abs(mu_all)));
            end
        else
            bif_types{iFw} = 'none';
        end

        prev_mu_all = mu_all;
    end

    bif_data{iOm} = struct('Om', Om, 'Fw_vec', Fw_vec, ...
        'A_resp', A_resp, 'maxMu_arr', maxMu_arr, ...
        'stable_arr', stable_arr, 'bif_types', {bif_types}, ...
        'mu_all_store', {mu_all_store});

    % Collect dominant multiplier trajectory for Argand diagram
    mu_traj = nan(N_Fw, 1);
    for iFw = 1:N_Fw
        if ~isempty(mu_all_store{iFw})
            [~, idx] = max(abs(mu_all_store{iFw}));
            mu_traj(iFw) = mu_all_store{iFw}(idx);
        end
    end
    mu_trajectories{iOm} = mu_traj;
end

%% ==================== Figure 1: Annotated Force Sweep ====================
figure('Color','w','Position',[50 50 1400 420]);
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for iOm = 1:3
    d = bif_data{iOm};
    nexttile;
    yyaxis left;
    ok_pts = isfinite(d.A_resp) & d.A_resp > 0;
    scatter(d.Fw_vec(d.stable_arr & ok_pts), d.A_resp(d.stable_arr & ok_pts), ...
        25, 'b', 'filled', 'DisplayName', 'Stable');
    hold on;
    scatter(d.Fw_vec(~d.stable_arr & ok_pts), d.A_resp(~d.stable_arr & ok_pts), ...
        25, 'r', 'DisplayName', 'Unstable');
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('F_w', 'FontName', fontName, 'FontSize', fsLab);
    ylabel('|x_1|_{amp}', 'FontName', fontName, 'FontSize', fsLab);
    grid on; box on;
    title(sprintf('\\Omega = %.2f', d.Om), 'FontName', fontName, 'FontSize', fsTit);

    % Annotate bifurcation types
    for iFw = 2:length(d.Fw_vec)
        if ok_pts(iFw) && ~strcmp(d.bif_types{iFw}, 'none')
            bif = d.bif_types{iFw};
            switch bif
                case 'fold'
                    label = 'Fold';
                    color = [1.0 0.2 0.2];
                case 'flip'
                    label = 'Flip';
                    color = [0.8 0.0 0.8];
                case 'neimark_sacker'
                    label = 'N-S';
                    color = [0.0 0.6 0.0];
                otherwise
                    label = bif;
                    color = [0.5 0.5 0.5];
            end
            text(d.Fw_vec(iFw), d.A_resp(iFw)*1.3, label, ...
                'FontName', fontName, 'FontSize', 8, 'Color', color, ...
                'FontWeight', 'bold');
        end
    end

    yyaxis right;
    plot(d.Fw_vec(ok_pts), d.maxMu_arr(ok_pts), 'k.-', 'MarkerSize', 6, 'LineWidth', 0.8);
    yline(tol_stable, 'k--', 'LineWidth', 1.0);
    yline(1.0, ':', 'Color', [0.4 0.4 0.4]);
    ylabel('max|\mu|', 'FontName', fontName, 'FontSize', fsLab);
    ylim([0.5, 2.5]);
    legend('Location', 'best', 'FontName', fontName, 'FontSize', 7);
end

sgtitle('Force Sweep: Bifurcation Type Classification (Fold / Flip / Neimark-Sacker)', ...
    'FontName', fontName, 'FontSize', 14);

% Save figure
out_dir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
exportgraphics(gcf, fullfile(out_dir, 'bifurcation_classification.pdf'), 'ContentType', 'vector');
fprintf('\n  -> bifurcation_classification.pdf saved\n');

%% ==================== Figure 2: Argand Diagram ====================
% Plot dominant multiplier trajectory in complex plane
figure('Color','w','Position',[60 60 1100 400]);
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for iOm = 1:3
    d = bif_data{iOm};
    nexttile;
    hold on; box on; grid on; axis equal;

    % Draw unit circle
    theta_c = linspace(0, 2*pi, 200);
    plot(cos(theta_c), sin(theta_c), 'k-', 'LineWidth', 1.5);

    % Plot dominant multiplier trajectory
    ok_pts = isfinite(d.A_resp) & d.A_resp > 0;
    for iFw = 1:length(d.Fw_vec)
        if ~ok_pts(iFw), continue; end
        mus = d.mu_all_store{iFw};
        if isempty(mus), continue; end

        % Plot all multipliers for this point
        for j = 1:length(mus)
            if d.stable_arr(iFw)
                plot(real(mus(j)), imag(mus(j)), 'b.', 'MarkerSize', 6);
            else
                plot(real(mus(j)), imag(mus(j)), 'r.', 'MarkerSize', 6);
            end
        end

        % Highlight dominant one
        [~, idx] = max(abs(mus));
        if d.stable_arr(iFw)
            plot(real(mus(idx)), imag(mus(idx)), 'bo', 'MarkerSize', 8, 'LineWidth', 1.0);
        else
            plot(real(mus(idx)), imag(mus(idx)), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
        end
    end

    % Mark special points
    plot(1, 0, 'k+', 'MarkerSize', 15, 'LineWidth', 2);
    plot(-1, 0, 'k+', 'MarkerSize', 15, 'LineWidth', 2);

    xlabel('Re(\mu)', 'FontName', fontName, 'FontSize', fsLab);
    ylabel('Im(\mu)', 'FontName', fontName, 'FontSize', fsLab);
    title(sprintf('\\Omega = %.2f', d.Om), 'FontName', fontName, 'FontSize', fsTit);
    xlim([-2.5, 2.5]); ylim([-2.0, 2.0]);
end

sgtitle('Floquet Multiplier Trajectories in Complex Plane (Argand Diagram)', ...
    'FontName', fontName, 'FontSize', 14);

exportgraphics(gcf, fullfile(out_dir, 'argand_multipliers.pdf'), 'ContentType', 'vector');
fprintf('  -> argand_multipliers.pdf saved\n');

%% ==================== Summary ====================
fprintf('\n========================================\n');
fprintf('  Bifurcation Classification Summary\n');
fprintf('========================================\n');
for iOm = 1:3
    d = bif_data{iOm};
    fprintf('\nOmega = %.2f:\n', d.Om);
    types_found = {};
    for iFw = 1:length(d.Fw_vec)
        if ~strcmp(d.bif_types{iFw}, 'none')
            types_found{end+1} = sprintf('Fw=%.4f: %s', d.Fw_vec(iFw), d.bif_types{iFw});
        end
    end
    if isempty(types_found)
        fprintf('  No bifurcation detected in the sweep range.\n');
    else
        for k = 1:length(types_found)
            fprintf('  %s\n', types_found{k});
        end
    end
end

%% ==================== Helper: Get Full Multiplier Set ====================
function mu_all = get_full_multipliers_fast(x_coeff, sysP, Omega, Nt)
% Compute all 6 Floquet multipliers (complex eigenvalues of monodromy matrix)
% Uses the same RK4 integrator as compute_floquet_fast but returns full set.
    x1c = x_coeff(1:5);
    x2c = x_coeff(6:10);
    qc  = x_coeff(11:15);

    Tp = 2*pi/Omega;
    dt = Tp / Nt;

    A_func = @(t) build_A_matrix_local(t, sysP, Omega, x1c, x2c, qc);

    Phi = eye(6);
    t = 0;
    for i = 1:Nt
        A = A_func(t);
        k1 = A * Phi;
        A_mid = A_func(t+0.5*dt);
        k2 = A_mid * (Phi + 0.5*dt*k1);
        k3 = A_mid * (Phi + 0.5*dt*k2);
        A_end = A_func(t+dt);
        k4 = A_end * (Phi + dt*k3);
        Phi = Phi + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
        t = t + dt;
    end
    mu_all = eig(Phi);
end

function A = build_A_matrix_local(t, sysP, Omega, x1c, x2c, qc)
    w = Omega;
    ct=cos(w*t); st=sin(w*t); c3t=cos(3*w*t); s3t=sin(3*w*t);

    x1 = x1c(1)+x1c(2)*ct+x1c(3)*st+x1c(4)*c3t+x1c(5)*s3t;
    x2 = x2c(1)+x2c(2)*ct+x2c(3)*st+x2c(4)*c3t+x2c(5)*s3t;

    be1=sysP(1); be2=sysP(2); mu=sysP(3);
    al1=sysP(4); ga1=sysP(5); ze=sysP(6);
    lam=sysP(7); kap_e=sysP(8); kap_c=sysP(9); sigma=sysP(10); ga2=sysP(11);

    theta = sqrt(max(lam,0));
    dx = x1-x2;

    df12 = (be1+al1) + 3*ga1*dx^2;
    df2g_x = be2 + 3*ga2*x2^2;
    df2g_v = 2*mu*ze;

    A = zeros(6);
    A(1,2) = 1;
    A(2,1) = -df12;       A(2,3) = df12;                        A(2,6) = theta;
    A(3,4) = 1;
    A(4,1) = df12/mu;     A(4,3) = (-df12 - df2g_x)/mu;         A(4,4) = -df2g_v/mu;  A(4,6) = -theta/mu;
    A(5,6) = 1;

    if abs(kap_e) > 1e-14
        A(6,2) = -theta/kap_e;   A(6,4) = theta/kap_e;
        A(6,5) = -kap_c/kap_e;   A(6,6) = -sigma/kap_e;
    else
        s = max(abs(sigma), 1e-12);
        A(6,2) = -50*theta/s;    A(6,4) = 50*theta/s;
        A(6,5) = -50*kap_c/s;    A(6,6) = -50;
    end
end
