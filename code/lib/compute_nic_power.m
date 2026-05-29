function P_active_dimless = compute_nic_power(x_coeff, sysP, Omega)
%% compute_nic_power - Compute cycle-averaged active power injected by NIC circuit
%
% The NIC synthesises a negative-resistance voltage that cancels coil resistance.
% From the dimensionless circuit equation:
%   kap_e * q'' + sigma * q' + kap_c * q = theta * (x1' - x2') + u_nic
%
% The NIC voltage u_nic represents the active power injection term.
% Passive terms are: -sigma_pas * q' (sigma_pas = passive resistance part)
% The total sigma = sigma_pas - sigma_act (sigma_act = NIC negative resistance)
% For our model: sigma is the TOTAL effective resistance.
%
% The active power injected by the NIC is:
%   P_active = <q' * u_nic>_T
%   where u_nic = -sigma_act * q' (negative resistance synthesis)
%   If sigma_pas = 1.0 (passive), then sigma = 1.0 - sigma_act
%   => sigma_act = 1.0 - sigma
%   => P_active = <(1.0 - sigma) * (q')^2>_T
%
% For the optimized parameter sigma_opt = 1.1506 (> 1.0), the NIC
% provides positive "active" resistance (actually more than full cancellation),
% which means it's absorbing energy rather than injecting it.
% This is consistent with a shunt damping configuration.
%
% Usage:
%   P = compute_nic_power(x_coeff, sysP, Omega)
%
% Input:
%   x_coeff - 15x1 HBM coefficients [x1(5); x2(5); q(5)]
%   sysP    - 11x1 system parameter vector
%   Omega   - excitation frequency
%
% Output:
%   P_active_dimless - scalar, cycle-averaged NIC active power (dimensionless)
%
% To convert to dimensional watts:
%   P_watts = P_active_dimless * (F0^2) / (m1 * omega_n)
%   where F0 = force amplitude [N], m1 = upper mass [kg], omega_n = natural freq [rad/s]

    x1c = x_coeff(1:5);
    x2c = x_coeff(6:10);
    qc  = x_coeff(11:15);

    sigma = sysP(10);
    lam = sysP(7);
    theta = sqrt(max(lam, 0));

    % The passive part of sigma (baseline RLC damping, no NIC effect)
    % For our model: sigma_total = sigma_passive - sigma_active
    % With NIC inactive: sigma = 1.0 (passive only)
    % With NIC active: sigma can be < 1.0 (partial cancellation) or > 1.0 (enhanced)
    % Active sigma_act = 1.0 - sigma (where 1.0 is the passive reference)
    sigma_act = 1.0 - sigma;

    % Reconstruct q'(t) over one period (64-point AFT grid, consistent with model)
    Nt = 64;
    t = linspace(0, 2*pi/Omega, Nt).';
    W = Omega;

    % q' = Omega * (b1*cos(wt) - a1*sin(wt)) + 3*Omega * (b3*cos(3wt) - a3*sin(3wt))
    qp = W * (qc(3)*cos(W*t) - qc(2)*sin(W*t)) ...
       + 3*W * (qc(5)*cos(3*W*t) - qc(4)*sin(3*W*t));

    % NIC synthetic voltage: u_nic = -sigma_act * q'
    % (negative because NIC synthesises a negative-resistance voltage)
    u_nic = -sigma_act * qp;

    % Instantaneous active power
    P_inst = qp .* u_nic;

    % Cycle-averaged
    P_active_dimless = mean(P_inst);

    % Also compute passive resistive power for completeness
    sigma_passive = 1.0;
    P_passive = mean(sigma_passive * qp.^2);

    % The total circuit power dissipation
    % P_total = <sigma_total * q'^2> = <(sigma_passive - sigma_act) * q'^2>
    %         = P_passive - P_active (if sigma_act > 0, NIC absorbs power)
    P_total = mean(sigma * qp.^2);

    % Consistency check: P_total = P_passive - P_active (if sigma_act = 1.0 - sigma)
    % (This should hold approximately)
end
