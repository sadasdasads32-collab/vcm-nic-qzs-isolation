function [bif_type, crossing_detail] = classify_bifurcation(mu_all_prev, mu_all_curr, tol)
%% classify_bifurcation - Classify bifurcation type from Floquet multipliers
%
% Detects how the dominant multiplier(s) cross the unit circle between
% two adjacent continuation steps, classifying into:
%   'fold'           - Real multiplier exits at (+1, 0)
%   'flip'            - Real multiplier exits at (-1, 0)
%   'neimark_sacker'  - Complex conjugate pair exits unit circle
%   'multiple'        - More than one mode detected simultaneously
%   'none'            - No stability change detected, or threshold not crossed
%
% Usage:
%   [bif_type, detail] = classify_bifurcation(mu_prev, mu_curr)
%   [bif_type, detail] = classify_bifurcation(mu_prev, mu_curr, tol)
%
% Input:
%   mu_all_prev  - Nx1 complex Floquet multipliers at previous step
%   mu_all_curr  - Nx1 complex Floquet multipliers at current step
%   tol          - Stability threshold (default 1.002)
%
% Output:
%   bif_type        - 'none' | 'fold' | 'flip' | 'neimark_sacker' | 'multiple'
%   crossing_detail - struct with fields:
%       .max_mu_prev, .max_mu_curr  - max |mu| at each step
%       .stable_prev, .stable_curr  - boolean stability flags
%       .fold_detected, .flip_detected, .ns_detected - booleans
%       .fold_mult, .flip_mult, .ns_mults - the relevant multipliers
%       .distance_change - change in distance to unit circle

    if nargin < 3 || isempty(tol)
        tol = 1.002;
    end

    crossing_detail = struct();
    crossing_detail.max_mu_prev = max(abs(mu_all_prev));
    crossing_detail.max_mu_curr = max(abs(mu_all_curr));
    crossing_detail.stable_prev = crossing_detail.max_mu_prev < tol;
    crossing_detail.stable_curr = crossing_detail.max_mu_curr < tol;

    % No stability change
    if crossing_detail.stable_prev == crossing_detail.stable_curr
        bif_type = 'none';
        crossing_detail.fold_detected = false;
        crossing_detail.flip_detected = false;
        crossing_detail.ns_detected = false;
        crossing_detail.fold_mult = [];
        crossing_detail.flip_mult = [];
        crossing_detail.ns_mults = [];
        crossing_detail.distance_change = 0;
        return;
    end

    %% Analyse which multiplier(s) crossed the unit circle

    % Find closest-to-real multipliers and complex pairs
    N = length(mu_all_curr);
    abs_mu = abs(mu_all_curr);

    % --- Fold detection: real multiplier crossing at (+1, 0) ---
    fold_detected = false;
    fold_mult = [];
    for i = 1:N
        mu = mu_all_curr(i);
        % Real multiplier (imaginary part very small) with magnitude > tol
        if abs(imag(mu)) < 0.05 * abs(real(mu)) && abs_mu(i) >= tol
            if real(mu) > 0  % crossing at (+1, 0)
                fold_detected = true;
                fold_mult = mu;
                break;
            end
        end
    end

    % --- Flip detection: real multiplier crossing at (-1, 0) ---
    flip_detected = false;
    flip_mult = [];
    for i = 1:N
        mu = mu_all_curr(i);
        if abs(imag(mu)) < 0.05 * abs(real(mu)) && abs_mu(i) >= tol
            if real(mu) < 0  % crossing at (-1, 0)
                flip_detected = true;
                flip_mult = mu;
                break;
            end
        end
    end

    % --- Neimark-Sacker detection: complex pair outside unit circle ---
    ns_detected = false;
    ns_mults = [];
    for i = 1:N
        mu = mu_all_curr(i);
        if abs(imag(mu)) >= 0.05 * abs(real(mu)) && abs_mu(i) >= tol
            % Find its conjugate pair
            for j = (i+1):N
                if abs(conj(mu) - mu_all_curr(j)) < 1e-8
                    ns_detected = true;
                    ns_mults = [mu; mu_all_curr(j)];
                    break;
                end
            end
            if ns_detected, break; end
        end
    end

    % --- Classify ---
    detected = [fold_detected, flip_detected, ns_detected];
    if sum(detected) == 0
        bif_type = 'none';
    elseif sum(detected) == 1
        if fold_detected
            bif_type = 'fold';
        elseif flip_detected
            bif_type = 'flip';
        else
            bif_type = 'neimark_sacker';
        end
    else
        bif_type = 'multiple';
    end

    crossing_detail.fold_detected = fold_detected;
    crossing_detail.flip_detected = flip_detected;
    crossing_detail.ns_detected = ns_detected;
    crossing_detail.fold_mult = fold_mult;
    crossing_detail.flip_mult = flip_mult;
    crossing_detail.ns_mults = ns_mults;
    crossing_detail.distance_change = crossing_detail.max_mu_curr - crossing_detail.max_mu_prev;
end
