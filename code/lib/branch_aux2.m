function [y, Jac] = branch_aux2(x, sysP)
%% branch_aux2
% 16维增广残差：
%
%   y = [ R(x); g(x) ]
%
% 其中：
%   R(x) : 15维 HBM 残差
%   g(x) : 1维加权切向伪弧长约束
%
% 伪弧长约束：
%
%   g = tangent^T * W * (x - xc) - arc = 0
%
% 其中：
%   xc      : 当前已知分支点
%   tangent : 当前分支切向量，已按加权范数归一化
%   arc     : 当前弧长步长
%   W       : diag([1,...,1, wp_arc^2])
%
% 注意：
%   x 是 16x1，其中 1:15 是 HBM 系数，第16维是延拓参数 Omega 或 F。
%
% 依赖 global:
%   tracking_file_name
%   xc
%   arc
%   tangent
%   wp_arc

    global tracking_file_name xc arc tangent wp_arc

    x = x(:);
    xc = xc(:);
    tangent = tangent(:);

    if numel(x) ~= 16
        error('branch_aux2: x must be 16x1.');
    end

    if isempty(xc) || numel(xc) ~= 16
        error('branch_aux2: xc is empty or not 16x1.');
    end

    if isempty(tangent) || numel(tangent) ~= 16
        error('branch_aux2: tangent is empty or not 16x1.');
    end

    if isempty(wp_arc)
        wp_arc = 15;
    end

    if isempty(arc) || ~isfinite(arc) || arc <= 0
        error('branch_aux2: arc is empty or invalid.');
    end

    %% ---------------------------------------------------------
    % 1) HBM 主残差及其对状态变量的 Jacobian
    %% ---------------------------------------------------------
    % tracking_file_name 通常为 'nondim_temp2'
    % 对 15维完整模型：
    %   R  : 15x1
    %   Jx : 15x15，对前15个 HBM 系数求导
    [R, Jx] = feval(tracking_file_name, x, sysP);

    if numel(R) ~= 15
        error('branch_aux2: tracking residual must be 15x1.');
    end

    if nargout > 1
        if isempty(Jx) || any(size(Jx) ~= [15,15])
            error('branch_aux2: tracking Jacobian must be 15x15.');
        end

        %% -----------------------------------------------------
        % 2) 对第16维延拓参数做有限差分，补齐 15x16 Jacobian
        %% -----------------------------------------------------
        h = 1e-7 * max(1, abs(x(16)));

        xt = x;
        xt(16) = xt(16) + h;

        Rt = feval(tracking_file_name, xt, sysP);

        dRdp = (Rt - R) / h;

        Jfull = [Jx, dRdp];
    end

    %% ---------------------------------------------------------
    % 3) 加权切向伪弧长约束
    %% ---------------------------------------------------------
    % W = diag([1,...,1, wp_arc^2])
    % 为了避免显式构造大矩阵，这里直接写成加权内积。

    dx = x - xc;

    g = tangent(1:15).' * dx(1:15) ...
        + (wp_arc^2) * tangent(16) * dx(16) ...
        - arc;

    y = [R; g];

    %% ---------------------------------------------------------
    % 4) 增广 Jacobian
    %% ---------------------------------------------------------
    if nargout > 1
        Jg = [tangent(1:15).', (wp_arc^2)*tangent(16)];

        Jac = [Jfull;
               Jg];
    end
end