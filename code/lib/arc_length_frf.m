function [Om, TF_dB, x_res] = arc_length_frf(sysP, Omega_Start, varargin)
%% arc_length_frf - Reusable arc-length FRF sweep with AFT cubic projection
%
% Computes the force transmissibility FRF using arc-length continuation
% of the HBM residuum for the electromechanical system.
%
% Usage:
%   [Om, TF_dB, x_res] = arc_length_frf(sysP, Omega_Start)
%   [Om, TF_dB, x_res] = arc_length_frf(sysP, Omega_Start, 'Fw', 0.005)
%   [Om, TF_dB, x_res] = arc_length_frf(sysP, Omega_Start, 'Step', -0.01, 'Steps', 3000)
%   [Om, TF_dB, x_res] = arc_length_frf(..., 'Plot', true)
%
% Input:
%   sysP        - 11x1 system parameter vector [be1,be2,mu,al1,ga1,ze1,lam,kap_e,kap_c,sigma,ga2]
%   Omega_Start - Starting frequency for arc-length continuation
%
% Optional PV pairs:
%   'Fw', val      - Excitation amplitude (default: global Fw or 0.005)
%   'Step', val    - Arc-length step size (default: -0.01, negative = downward)
%   'Steps', val   - Number of continuation steps (default: 3000)
%   'Plot', tf     - Whether to plot the FRF (default: false)
%   'DisplayName', str - Legend entry for plot (default: 'FRF')
%
% Output:
%   Om    - Frequency vector [Nx1]
%   TF_dB - Force transmissibility in dB [Nx1], 20*log10(|ft|/Fw)
%   x_res - Full branch follow result [16xN] (last row = Omega)

    %% Parse inputs
    p = inputParser;
    p.addRequired('sysP', @(x) isnumeric(x) && numel(x)==11);
    p.addRequired('Omega_Start', @(x) isscalar(x) && x>0);
    p.addParameter('Fw', [], @(x) isscalar(x) && x>0);
    p.addParameter('Step', -0.01, @(x) isscalar(x));
    p.addParameter('Steps', 3000, @(x) isscalar(x) && x>=100);
    p.addParameter('Plot', false, @islogical);
    p.addParameter('DisplayName', 'FRF', @ischar);
    p.parse(sysP, Omega_Start, varargin{:});

    opts = p.Results;
    Omega_Step = opts.Step;
    N_steps    = opts.Steps;

    % Determine Fw
    global Fw
    if ~isempty(opts.Fw)
        Fw_use = opts.Fw;
    elseif ~isempty(Fw) && isfinite(Fw)
        Fw_use = Fw;
    else
        Fw_use = 0.005;
    end

    %% Solve first two points via Newton
    Omega_Next = Omega_Start + Omega_Step;

    % Point 0: at Omega_Start
    y_init = zeros(15, 1);
    y_init(end+1) = Omega_Start;
    global FixedOmega
    FixedOmega_backup = FixedOmega;
    FixedOmega = [];
    Fw_backup = Fw;
    Fw = Fw_use;

    [x0_full, ok0] = newton('nondim_temp2', y_init, sysP);
    if ~ok0
        warning('arc_length_frf: Newton failed at Omega_Start=%.4f', Omega_Start);
    end
    x0 = x0_full(1:15);

    % Point 1: at Omega_Next
    y_init2 = [x0; Omega_Next];
    [x1_full, ok1] = newton('nondim_temp2', y_init2, sysP);
    if ~ok1
        warning('arc_length_frf: Newton failed at Omega_Next=%.4f', Omega_Next);
    end
    x1 = x1_full(1:15);

    %% Arc-length continuation
    [x_res, ~] = branch_follow2('nondim_temp2', N_steps, ...
                                Omega_Start, Omega_Next, x0, x1, sysP);

    %% Compute force transmissibility
    Om = x_res(16,:).';
    be2 = sysP(2);
    mu  = sysP(3);
    ze2 = sysP(6);
    ga2 = sysP(11);

    % Extract x2 coefficients [x20, a21, b21, a23, b23] (rows 6:10)
    x2 = x_res(6:10,:).';

    % Compute x2_dot in harmonic domain
    W = Om;
    x2_dot = zeros(size(x2));
    x2_dot(:,1) = 0;                    % dc component
    x2_dot(:,2) = W .* x2(:,3);        % cos1 -> omega*sin1
    x2_dot(:,3) = -W .* x2(:,2);       % sin1 -> -omega*cos1
    x2_dot(:,4) = 3*W .* x2(:,5);      % cos3 -> 3*omega*sin3
    x2_dot(:,5) = -3*W .* x2(:,4);     % sin3 -> -3*omega*cos3

    % AFT for cubic nonlinearity
    x2_cub = cubic_proj_013(x2);

    % Transmitted force
    ft = be2*x2 + ga2*x2_cub + 2*mu*ze2*x2_dot;
    ft1 = hypot(ft(:,2), ft(:,3));       % fundamental amplitude
    ft3 = hypot(ft(:,4), ft(:,5));       % 3rd harmonic amplitude
    ft_amp = hypot(ft1, ft3);             % total amplitude

    % Convert to dB
    TF = ft_amp ./ Fw_use;
    TF_dB = 20*log10(max(TF, 1e-300));

    % Filter invalid points
    valid = isfinite(Om) & isfinite(TF_dB) & (Om > 0);
    Om = Om(valid);
    TF_dB = TF_dB(valid);

    %% Plot (optional)
    if opts.Plot
        figure('Color','w', 'Position',[150 150 700 500]);
        ax = gca; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        semilogx(Om, TF_dB, 'b-', 'LineWidth', 1.5, ...
                 'DisplayName', opts.DisplayName);
        yline(0, 'k--', '0 dB');
        xlabel('\Omega');
        ylabel('T_F (dB)');
        title(sprintf(['BG Model FRF: \\lambda=%.2f, \\kappa_e=%.2f, ' ...
              '\\kappa_c=%.2f, \\sigma=%.2f'], ...
              sysP(7), sysP(8), sysP(9), sysP(10)));
        legend('Location', 'best');
        xlim([0.1, Omega_Start]);
    end

    %% Restore globals
    FixedOmega = FixedOmega_backup;
    Fw = Fw_backup;
end

%% ==================== Local AFT Helpers ====================
function cubic = cubic_proj_013(U)
% U: N x 5 matrix, each row = [dc, cos1, sin1, cos3, sin3]
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
        c1 = cos(t); s1 = sin(t); c3 = cos(3*t); s3 = sin(3*t);
        dc = ones(N, 1);
        pT = [dc, c1, s1, c3, s3];
        Inv = [dc, 2*c1, 2*s1, 2*c3, 2*s3]';
        pTinv = (1/N) * Inv;
        pTinv(1,:) = (1/N) * dc';
    end
    T_mat = pT; T_inv = pTinv;
end
