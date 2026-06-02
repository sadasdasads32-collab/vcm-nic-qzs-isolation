function [Om, TF_dB, x_res] = arc_length_frf_robust(sysP, Omega_Start, varargin)
%% arc_length_frf_robust - Robust multi-segment arc-length FRF sweep
%
% Uses upward sweep (0.1 -> 10.0) with multi-segment restart, designed
% for softening Duffing systems where downward sweep misses resonance.
%
% Drop-in replacement for arc_length_frf with same output API.
%
% Usage:
%   [Om, TF_dB, x_res] = arc_length_frf_robust(sysP, 10.0)
%   [Om, TF_dB, x_res] = arc_length_frf_robust(sysP, 10.0, 'Fw', 0.005)
%   [Om, TF_dB, x_res] = arc_length_frf_robust(..., 'Budget', 3000)
%   [Om, TF_dB, x_res] = arc_length_frf_robust(..., 'Verbose', false)
%
% Optional PV pairs:
%   'Fw', val      - Excitation amplitude (default: global Fw or 0.005)
%   'Budget', val  - Max total continuation steps (default: 8000)
%   'Verbose', tf  - Print progress (default: true)
%
% Output:
%   Om    - Frequency vector [Nx1]
%   TF_dB - Force transmissibility in dB [Nx1]
%   x_res - Full branch follow result [16xN] (last row = Omega)

    %% Parse inputs
    p = inputParser;
    p.addRequired('sysP', @(x) isnumeric(x) && numel(x)==11);
    p.addRequired('Omega_Start', @(x) isscalar(x) && x>0);  % kept for API compat
    p.addParameter('Fw', [], @(x) isscalar(x) && x>0);
    p.addParameter('Budget', 8000, @(x) isscalar(x) && x>=200);
    p.addParameter('Verbose', true, @islogical);
    p.parse(sysP, Omega_Start, varargin{:});
    opts = p.Results;

    %% Determine Fw
    global Fw
    if ~isempty(opts.Fw)
        Fw_use = opts.Fw;
    elseif ~isempty(Fw) && isfinite(Fw)
        Fw_use = Fw;
    else
        Fw_use = 0.005;
    end

    %% Save and set globals
    global FixedOmega ParamMin ParamMax
    Fw_backup = Fw;
    FixedOmega_backup = FixedOmega;
    ParamMin_backup = ParamMin;
    ParamMax_backup = ParamMax;

    Fw = Fw_use;
    FixedOmega = [];

    %% Multi-segment upward FRF sweep
    try
        x_res = frf_upward_sweep(sysP, opts.Budget, opts.Verbose);
    catch ME
        Fw = Fw_backup;
        FixedOmega = FixedOmega_backup;
        ParamMin = ParamMin_backup;
        ParamMax = ParamMax_backup;
        rethrow(ME);
    end

    %% Restore globals
    Fw = Fw_backup;
    FixedOmega = FixedOmega_backup;
    ParamMin = ParamMin_backup;
    ParamMax = ParamMax_backup;

    %% Compute force transmissibility
    if isempty(x_res) || size(x_res,2) < 2
        Om = []; TF_dB = []; x_res = [];
        return;
    end

    Om = x_res(16,:).';
    be2 = sysP(2); mu = sysP(3); ze2 = sysP(6); ga2 = sysP(11);

    x2 = x_res(6:10,:).';
    W = Om;
    x2_dot = zeros(size(x2));
    x2_dot(:,1) = 0;
    x2_dot(:,2) = W .* x2(:,3);
    x2_dot(:,3) = -W .* x2(:,2);
    x2_dot(:,4) = 3*W .* x2(:,5);
    x2_dot(:,5) = -3*W .* x2(:,4);

    x2_cub = cubic_proj_013(x2);
    ft = be2*x2 + ga2*x2_cub + 2*mu*ze2*x2_dot;
    ft1 = hypot(ft(:,2), ft(:,3));
    ft3 = hypot(ft(:,4), ft(:,5));
    ft_amp = hypot(ft1, ft3);

    TF = ft_amp ./ Fw_use;
    TF_dB = 20*log10(max(TF, 1e-300));

    valid = isfinite(Om) & isfinite(TF_dB) & (Om > 0);
    Om = Om(valid);
    TF_dB = TF_dB(valid);

    if opts.Verbose
        fprintf('Robust FRF: Omega [%.4f, %.4f], %d valid points\n', ...
            min(Om), max(Om), length(Om));
    end
end

%% ==================== Internal: Upward Sweep Engine ====================
function x_all = frf_upward_sweep(sysP, budget, verbose)
% Multi-segment upward arc-length sweep from low to high frequency.
% Designed for softening systems: starts in linear regime and sweeps up.

    global ParamMin ParamMax

    Omega0    = 0.1;
    OmegaMax  = 10.0;
    chunk_steps = min(1200, budget);
    dOmega_min  = 1e-5;
    init_eps    = 1e-4;
    min_pts_ok  = 20;

    ParamMin = Omega0 - 0.05;
    ParamMax = OmegaMax;

    %% Point 0: Newton at Omega0 (near zero, linear regime)
    y0 = zeros(15,1);
    y0(2) = init_eps;
    [x0_full, ok0] = newton_safe('nondim_temp2', [y0; Omega0], sysP);
    if ~ok0
        if verbose
            warning('FRF_robust: Newton failed at Omega0=%.4f', Omega0);
        end
        x_all = [];
        return;
    end
    x0 = x0_full(1:15);

    %% Point 1: small step upward to set tangent direction
    dOm = 1e-2;
    ok1 = false;
    while dOm >= dOmega_min
        Omega1 = Omega0 + dOm;
        [x1_full, ok1] = newton_safe('nondim_temp2', [x0; Omega1], sysP);
        if ok1, break; end
        dOm = dOm * 0.5;
    end
    if ~ok1
        if verbose
            warning('FRF_robust: Cannot construct 2nd point at Omega0=%.4f', Omega0);
        end
        x_all = [];
        return;
    end
    x1 = x1_full(1:15);

    if verbose
        fprintf('FRF_robust: Omega0=%.4f, Omega1=%.4f (dOm=%.1e)\n', Omega0, Omega1, dOm);
    end

    %% Multi-segment continuation
    x_all = [];
    cur_x0 = x0; cur_x1 = x1;
    cur_Om0 = Omega0; cur_Om1 = Omega1;
    remain_budget = budget;

    while remain_budget > 0
        nsteps = min(chunk_steps, remain_budget);

        try
            [x_seg, ~] = branch_follow2('nondim_temp2', nsteps, ...
                cur_Om0, cur_Om1, cur_x0, cur_x1, sysP);
        catch ME
            if verbose
                fprintf('FRF_robust: branch_follow2 crashed: %s\n', ME.message);
            end
            break;
        end

        if isempty(x_seg) || size(x_seg,1) ~= 16 || size(x_seg,2) < 2
            break;
        end

        % Trim NaN/Inf
        Om_seg = x_seg(16,:);
        good = isfinite(Om_seg);
        if ~all(good)
            last_good = find(good, 1, 'last');
            if isempty(last_good) || last_good < 2, break; end
            x_seg = x_seg(:, 1:last_good);
        end

        % Concatenate (avoid duplicating first point)
        if isempty(x_all)
            x_all = x_seg;
        else
            x_all = [x_all, x_seg(:, 2:end)]; %#ok<AGROW>
        end

        % Check termination
        if x_all(16,end) >= OmegaMax
            if verbose
                fprintf('FRF_robust: Reached OmegaMax=%.3f\n', OmegaMax);
            end
            break;
        end

        % If this segment produced too few points, try smaller restart
        if size(x_seg,2) < min_pts_ok
            cur_Om0 = x_all(16,end);
            cur_x0  = x_all(1:15,end);
            dOm2 = 5e-3; ok2 = false;
            while dOm2 >= dOmega_min
                cur_Om1 = cur_Om0 + dOm2;
                [tmp, ok2] = newton_safe('nondim_temp2', [cur_x0; cur_Om1], sysP);
                if ok2, break; end
                dOm2 = dOm2 * 0.5;
            end
            if ~ok2, break; end
            cur_x1 = tmp(1:15);
        else
            cur_Om0 = x_all(16,end-1);
            cur_x0  = x_all(1:15,end-1);
            cur_Om1 = x_all(16,end);
            cur_x1  = x_all(1:15,end);
        end

        remain_budget = budget - size(x_all,2);
    end

    %% Trim to OmegaMax
    if ~isempty(x_all)
        Om_final = x_all(16,:);
        idx = find(Om_final <= OmegaMax);
        if ~isempty(idx)
            x_all = x_all(:, 1:idx(end));
        end
    end
end

%% ==================== Safe Newton Wrapper ====================
function [x, ok] = newton_safe(funname, x0, sysP)
    x = []; ok = false;
    try
        [x, ok] = newton(funname, x0, sysP);
    catch
        x = []; ok = false;
    end
end

%% ==================== AFT Cubic Projection ====================
function cubic = cubic_proj_013(U)
    [T_mat, T_inv] = get_AFT_matrices();
    X_time  = (T_mat * U.').';
    X3_time = X_time.^3;
    cubic   = (T_inv * X3_time.').';
end

function [T_mat, T_inv] = get_AFT_matrices()
    persistent pT pTinv
    if isempty(pT)
        N = 64;
        t = (0:N-1)'*(2*pi/N);
        c1=cos(t); s1=sin(t); c3=cos(3*t); s3=sin(3*t); dc=ones(N,1);
        pT = [dc, c1, s1, c3, s3];
        Inv = [dc, 2*c1, 2*s1, 2*c3, 2*s3]';
        pTinv = (1/N) * Inv;
        pTinv(1,:) = (1/N) * dc';
    end
    T_mat = pT; T_inv = pTinv;
end
