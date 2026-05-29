function [x_branch, info] = branch_follow2(tracking_fun, max_steps, mu0, mu1, x0, x1, sysP)
%% branch_follow2
% 15维完整模型的加权切向伪弧长延拓
%
% 输入：
%   tracking_fun : 残差函数名，例如 'nondim_temp2'
%   max_steps    : 最大延拓步数
%   mu0          : 第一个点的参数值，例如 Omega_Start
%   mu1          : 第二个点的参数值，例如 Omega_Next
%   x0           : 第一个点的 15维 HBM 系数
%   x1           : 第二个点的 15维 HBM 系数
%   sysP         : 系统参数
%
% 输出：
%   x_branch : 16 x N
%              前15行为 HBM 系数，第16行为延拓参数
%   info     : 结构体，记录延拓信息
%
% 说明：
%   本版本使用真正的切向伪弧长约束：
%
%       tangent^T W (x - xc) - arc = 0
%
%   而不是 norm(x-xc)-arc=0。
%
%   该形式更适合追踪折叠分支、多值响应。
%
% 依赖：
%   newton.m
%   branch_aux2.m
%   tracking_fun，例如 nondim_temp2.m
%
% 需要 global:
%   ParamMin, ParamMax 可选，用来限定参数范围。

    global tracking_file_name xc arc tangent wp_arc
    global ParamMin ParamMax

    %% ---------------------------------------------------------
    % 0) 基本检查
    %% ---------------------------------------------------------
    if numel(x0) ~= 15 || numel(x1) ~= 15
        error('branch_follow2: x0 and x1 must be 15x1.');
    end

    x0 = x0(:);
    x1 = x1(:);

    if nargin < 2 || isempty(max_steps)
        max_steps = 3000;
    end

    if isempty(ParamMin)
        ParamMin = min(mu0, mu1) - 0.05;
    end

    if isempty(ParamMax)
        ParamMax = max(mu0, mu1) + 0.05;
    end

    tracking_file_name = tracking_fun;

    %% ---------------------------------------------------------
    % 1) 用户可调参数
    %% ---------------------------------------------------------
    wp_arc = 15;          % 参数方向加权，越大表示 Omega 变化被放大
    arc_min = 1e-5;       % 最小弧长
    arc_max = 0.15;       % 最大弧长
    arc_grow = 1.15;      % 成功后放大步长
    arc_shrink = 0.5;     % 失败后缩小步长

    max_try = 10;         % 每一步最大尝试次数
    newton_tol_accept = 1e-6;

    print_every = 50;

    %% ---------------------------------------------------------
    % 2) 初始两个完整点
    %% ---------------------------------------------------------
    X0 = [x0; mu0];
    X1 = [x1; mu1];

    x_branch = [X0, X1];

    %% ---------------------------------------------------------
    % 3) 初始切向量
    %% ---------------------------------------------------------
    dX = X1 - X0;

    arc0 = weighted_norm_local(dX, wp_arc);

    if arc0 <= 0 || ~isfinite(arc0)
        error('branch_follow2: initial two points are invalid.');
    end

    tangent = dX / arc0;

    % 初始弧长
    arc = min(max(arc0, arc_min), arc_max);

    fprintf('Branch Follow Pseudo-Arc: Start. Arc=%.2e, wp=%g, Range=[%.6f, %.6f]\n', ...
            arc, wp_arc, ParamMin, ParamMax);

    %% ---------------------------------------------------------
    % 4) 主循环
    %% ---------------------------------------------------------
    stop_reason = 'max_steps';
    fail_count_total = 0;

    for step = 2:max_steps

        Xcur = x_branch(:,end);

        % 当前点作为伪弧长约束中的 xc
        xc = Xcur;

        success = false;
        arc_try = arc;

        %% -----------------------------------------------------
        % 4.1 当前步尝试
        %% -----------------------------------------------------
        for attempt = 1:max_try

            % Predictor
            Xpred = Xcur + arc_try * tangent;

            % 把当前尝试步长传给 branch_aux2
            arc = arc_try;

            try
                [Xnew, ok, Rn] = newton('branch_aux2', Xpred, sysP);
            catch ME
                ok = false;
                Rn = inf;
                Xnew = [];
                fprintf('   [Step %d Try %d] Newton crash: %s\n', ...
                        step, attempt, ME.message);
            end

            if ok && ~isempty(Xnew) && all(isfinite(Xnew)) && Rn < newton_tol_accept

                mu_new = Xnew(16);

                % 参数范围检查
                % 注意：这里不要求参数单调，只限制在总范围内。
                if mu_new < ParamMin || mu_new > ParamMax
                    % 如果已经到边界附近，则认为追踪结束。
                    if abs(mu_new - ParamMin) < 5*arc_try || abs(mu_new - ParamMax) < 5*arc_try
                        stop_reason = sprintf('reached param bound %.6f', mu_new);
                        fprintf('   End: reached param bound %.6f (Range=[%.6f, %.6f])\n', ...
                                mu_new, ParamMin, ParamMax);
                        info = build_info_local(stop_reason, step, fail_count_total, arc, wp_arc);
                        return;
                    end

                    % 否则缩小步长再试
                    arc_try = arc_try * arc_shrink;

                    if arc_try < arc_min
                        stop_reason = sprintf('param out of range %.6f', mu_new);
                        fprintf('   Stop: param out of range %.6f, arc too small.\n', mu_new);
                        info = build_info_local(stop_reason, step, fail_count_total, arc, wp_arc);
                        return;
                    end

                    continue;
                end

                success = true;
                break;
            else
                arc_try = arc_try * arc_shrink;

                if arc_try < arc_min
                    break;
                end
            end
        end

        %% -----------------------------------------------------
        % 4.2 当前步失败
        %% -----------------------------------------------------
        if ~success
            fail_count_total = fail_count_total + 1;

            fprintf('   [Fail] Step %d failed. arc=%.2e\n', step, arc_try);

            if arc_try < arc_min
                stop_reason = 'arc below arc_min';
                break;
            end

            % 不重置为纯参数方向，只缩小步长并继续
            arc = max(arc_try, arc_min);

            continue;
        end

        %% -----------------------------------------------------
        % 4.3 接受新点
        %% -----------------------------------------------------
        x_branch = [x_branch, Xnew]; %#ok<AGROW>

        %% -----------------------------------------------------
        % 4.4 更新切向量
        %% -----------------------------------------------------
        dXnew = Xnew - Xcur;

        nrm = weighted_norm_local(dXnew, wp_arc);

        if nrm <= 0 || ~isfinite(nrm)
            stop_reason = 'invalid tangent';
            break;
        end

        tangent_new = dXnew / nrm;

        % 保持切向方向连续
        if weighted_dot_local(tangent_new, tangent, wp_arc) < 0
            tangent_new = -tangent_new;
        end

        tangent = tangent_new;

        %% -----------------------------------------------------
        % 4.5 自适应步长
        %% -----------------------------------------------------
        if attempt <= 2
            arc = min(arc_try * arc_grow, arc_max);
        elseif attempt <= 5
            arc = arc_try;
        else
            arc = max(arc_try * arc_shrink, arc_min);
        end

        %% -----------------------------------------------------
        % 4.6 打印信息
        %% -----------------------------------------------------
        if step <= 5 || mod(step, print_every) == 0
            fprintf('   Step %5d: Param=%.6f | arc=%.2e | R=%.2e | try=%d\n', ...
                    step, Xnew(16), arc, Rn, attempt);
        end

        %% -----------------------------------------------------
        % 4.7 边界停止
        %% -----------------------------------------------------
        if Xnew(16) <= ParamMin || Xnew(16) >= ParamMax
            stop_reason = sprintf('reached param bound %.6f', Xnew(16));
            fprintf('   End: reached param bound %.6f (Range=[%.6f, %.6f])\n', ...
                    Xnew(16), ParamMin, ParamMax);
            break;
        end
    end

    %% ---------------------------------------------------------
    % 5) 输出信息
    %% ---------------------------------------------------------
    fprintf('Branch Follow Pseudo-Arc: Finished. Steps=%d, LastParam=%.6f\n', ...
            size(x_branch,2)-1, x_branch(16,end));

    info = build_info_local(stop_reason, size(x_branch,2)-1, fail_count_total, arc, wp_arc);
end

%% =========================================================
% 加权范数
%% =========================================================
function n = weighted_norm_local(x, wp)
    x = x(:);
    n = sqrt(sum(x(1:end-1).^2) + (wp^2)*x(end)^2);
end

%% =========================================================
% 加权内积
%% =========================================================
function v = weighted_dot_local(a, b, wp)
    a = a(:);
    b = b(:);
    v = a(1:end-1).'*b(1:end-1) + (wp^2)*a(end)*b(end);
end

%% =========================================================
% 输出 info
%% =========================================================
function info = build_info_local(stop_reason, steps, fail_count_total, arc, wp_arc)
    info = struct();
    info.stop_reason = stop_reason;
    info.steps = steps;
    info.fail_count_total = fail_count_total;
    info.final_arc = arc;
    info.wp_arc = wp_arc;
end