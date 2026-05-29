function [R_aug, Jac_aug] = branch_aux2N(x, sysP)
%% branch_aux2N
% N维通用伪弧长增广残差
%
% 用于 10维算子模型时：
%   x : 11x1 = [x1(5); x2(5); Omega_or_Fw]
%
% 增广残差：
%   R_aug = [R_main; g]
%
% 其中：
%   R_main : n维主残差，例如 nondim_temp2_op 返回 10x1
%   g      : 加权切向伪弧长约束
%
% 伪弧长约束：
%   g = tangent^T * W * (x - xc) - arc = 0
%
% global 变量由 branch_follow2N 设置：
%   tracking_file_name
%   xc
%   arc
%   tangent
%   wp_arc

    global tracking_file_name xc arc tangent wp_arc

    x = x(:);
    xc = xc(:);
    tangent = tangent(:);

    n1 = numel(x);
    n  = n1 - 1;

    if n <= 0
        error('branch_aux2N: invalid dimension.');
    end

    if isempty(xc) || numel(xc) ~= n1
        error('branch_aux2N: xc is empty or dimension mismatch.');
    end

    if isempty(tangent) || numel(tangent) ~= n1
        error('branch_aux2N: tangent is empty or dimension mismatch.');
    end

    if isempty(wp_arc)
        wp_arc = 15;
    end

    if isempty(arc) || ~isfinite(arc) || arc <= 0
        error('branch_aux2N: arc is empty or invalid.');
    end

    %% ---------------------------------------------------------
    % 1) 主残差与状态 Jacobian
    %% ---------------------------------------------------------
    % 对 nondim_temp2_op:
    %   R_main : 10x1
    %   J_state: 10x10，对前10个机械系数求导
    [R_main, J_state] = feval(tracking_file_name, x, sysP);

    R_main = R_main(:);

    if numel(R_main) ~= n
        error('branch_aux2N: main residual size mismatch. Expected %d, got %d.', ...
              n, numel(R_main));
    end

    %% ---------------------------------------------------------
    % 2) 主残差对最后一个参数的有限差分导数
    %% ---------------------------------------------------------
    if nargout > 1

        if isempty(J_state) || any(size(J_state) ~= [n,n])
            error('branch_aux2N: main Jacobian must be %dx%d.', n, n);
        end

        h = 1e-7 * max(1, abs(x(end)));

        xt = x;
        xt(end) = xt(end) + h;

        Rt = feval(tracking_file_name, xt, sysP);
        Rt = Rt(:);

        dRdp = (Rt - R_main) / h;

        J_full = [J_state, dRdp];
    end

    %% ---------------------------------------------------------
    % 3) 加权切向伪弧长约束
    %% ---------------------------------------------------------
    dx = x - xc;

    g = tangent(1:n).' * dx(1:n) ...
        + (wp_arc^2) * tangent(n1) * dx(n1) ...
        - arc;

    R_aug = [R_main; g];

    %% ---------------------------------------------------------
    % 4) 增广 Jacobian
    %% ---------------------------------------------------------
    if nargout > 1
        Jg = [tangent(1:n).', (wp_arc^2)*tangent(n1)];

        Jac_aug = [J_full;
                   Jg];
    end
end