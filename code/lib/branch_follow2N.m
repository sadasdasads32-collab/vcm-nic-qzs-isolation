function [x_branch, info] = branch_follow2N(fname, nsteps, mu0, mu1, x0, x1, sysP)
%% branch_follow2N
% N维通用加权切向伪弧长延拓
%
% 主要用于 10维复算子模型：
%   fname = 'nondim_temp2_op'
%   x0, x1 为 10x1 状态向量
%   输出 x_branch 为 11xN
%
% 也可用于其他 n维状态 + 1个延拓参数的系统。
%
% 输入：
%   fname  : 主残差函数名，例如 'nondim_temp2_op'
%   nsteps : 最大延拓步数
%   mu0    : 第一个点参数值，例如 Omega_Start
%   mu1    : 第二个点参数值，例如 Omega_Next
%   x0     : 第一个点状态，n x 1 或 n+1 x 1
%   x1     : 第二个点状态，n x 1 或 n+1 x 1
%   sysP   : 系统参数
%
% 输出：
%   x_branch : (n+1) x N
%   info     : 结构体
%
% 特点：
%   1) 使用切向伪弧长约束，不使用球面约束；
%   2) 不强制参数单调；
%   3) 失败时只缩小步长，不重置为纯参数方向；
%   4) 适合追踪折叠分支和一频多值响应。

    global tracking_file_name xc arc tangent wp_arc
    global ParamMin ParamMax

    tracking_file_name = fname;

    %% ---------------------------------------------------------
    % 0) 输入整理
    %% ---------------------------------------------------------
    x0 = x0(:);
    x1 = x1(:);

    if numel(x0) ~= numel(x1)
        error('branch_follow2N: x0 and x1 must have the same length.');
    end

    % 如果传入的是状态向量 n x 1，则补参数；
    % 如果传入的是完整向量 n+1 x 1，则直接使用。
    if abs(x0(end) - mu0) < 1e-12 && abs(x1(end) - mu1) < 1e-12
        X0 = x0;
        X1 = x1;
    else
        X0 = [x0; mu0];
        X1 = [x1; mu1];
    end

    n1 = numel(X0);
    n  = n1 - 1;

    if n <= 0
        error('branch_follow2N: invalid dimension.');
    end

    if numel(X1) ~= n1
        error('branch_follow2N: X0/X1 dimension mismatch.');
    end

    %% ---------------------------------------------------------
    % 1) 参数范围
    %% ---------------------------------------------------------
    if isempty(ParamMin) || ~isfinite(ParamMin)
        ParamMin = min(mu0,mu1) - 0.05;
    end

    if isempty(ParamMax) || ~isfinite(ParamMax)
        ParamMax = max(mu0,mu1) + 0.05;
    end

    %% ---------------------------------------------------------
    % 2) 用户可调参数
    %% ---------------------------------------------------------
    wp_arc = 15;          % 参数方向加权
    arc_min = 1e-5;       % 最小弧长
    arc_max = 0.15;       % 最大弧长
    arc_grow = 1.15;      % 成功后步长放大
    arc_shrink = 0.5;     % 失败后步长缩小

    max_try = 10;         % 每步最大尝试次数
    newton_tol_accept = 1e-6;

    print_every = 50;

    %% ---------------------------------------------------------
    % 3) 初始化分支
    %% ---------------------------------------------------------
    x_branch = [X0, X1];

    dX = X1 - X0;
    arc0 = weighted_normN_local(dX, wp_arc);

    if arc0 <= 0 || ~isfinite(arc0)
        error('branch_follow2N: initial two points are invalid.');
    end

    tangent = dX / arc0;

    arc = min(max(arc0, arc_min), arc_max);

    fprintf('Branch Follow N Pseudo-Arc: Start. n=%d, Arc=%.2e, wp=%g, Range=[%.6f, %.6f]\n', ...
            n, arc, wp_arc, ParamMin, ParamMax);

    %% ---------------------------------------------------------
    % 4) 主循环
    %% ---------------------------------------------------------
    stop_reason = 'max_steps';
    fail_count_total = 0;

    for step = 2:nsteps

        Xcur = x_branch(:,end);

        % 当前点作为伪弧长约束基点
        xc = Xcur;

        success = false;
        arc_try = arc;

        %% -----------------------------------------------------
        % 4.1 当前步尝试
        %% -----------------------------------------------------
        for attempt = 1:max_try

            % Predictor
            Xpred = Xcur + arc_try * tangent;

            % 当前尝试弧长传给 branch_aux2N
            arc = arc_try;

            try
                [Xnew, ok, Rn] = newton('branch_aux2N', Xpred, sysP);
            catch ME
                ok = false;
                Rn = inf;
                Xnew = [];
                fprintf('   [Step %d Try %d] Newton crash: %s\n', ...
                        step, attempt, ME.message);
            end

            if ok && ~isempty(Xnew) && all(isfinite(Xnew)) && Rn < newton_tol_accept

                mu_new = Xnew(end);

                % 总参数范围限制，不强制参数单调
                if mu_new < ParamMin || mu_new > ParamMax

                    if abs(mu_new - ParamMin) < 5*arc_try || abs(mu_new - ParamMax) < 5*arc_try
                        stop_reason = sprintf('reached param bound %.6f', mu_new);
                        fprintf('   End: reached param bound %.6f (Range=[%.6f, %.6f])\n', ...
                                mu_new, ParamMin, ParamMax);
                        info = build_infoN_local(stop_reason, step, fail_count_total, arc, wp_arc, n);
                        return;
                    end

                    arc_try = arc_try * arc_shrink;

                    if arc_try < arc_min
                        stop_reason = sprintf('param out of range %.6f', mu_new);
                        fprintf('   Stop: param out of range %.6f, arc too small.\n', mu_new);
                        info = build_infoN_local(stop_reason, step, fail_count_total, arc, wp_arc, n);
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

            % 不重置为纯参数方向，只缩小步长继续
            arc = max(arc_try, arc_min);

            continue;
        end

        %% -----------------------------------------------------
        % 4.3 接受新点
        %% -----------------------------------------------------
        x_branch = [x_branch, Xnew]; %#ok<AGROW>

        %% -----------------------------------------------------
        % 4.4 更新切向
        %% -----------------------------------------------------
        dXnew = Xnew - Xcur;

        nrm = weighted_normN_local(dXnew, wp_arc);

        if nrm <= 0 || ~isfinite(nrm)
            stop_reason = 'invalid tangent';
            break;
        end

        tangent_new = dXnew / nrm;

        % 保持方向连续
        if weighted_dotN_local(tangent_new, tangent, wp_arc) < 0
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
                    step, Xnew(end), arc, Rn, attempt);
        end

        %% -----------------------------------------------------
        % 4.7 边界停止
        %% -----------------------------------------------------
        if Xnew(end) <= ParamMin || Xnew(end) >= ParamMax
            stop_reason = sprintf('reached param bound %.6f', Xnew(end));
            fprintf('   End: reached param bound %.6f (Range=[%.6f, %.6f])\n', ...
                    Xnew(end), ParamMin, ParamMax);
            break;
        end
    end

    %% ---------------------------------------------------------
    % 5) 输出 info
    %% ---------------------------------------------------------
    fprintf('Branch Follow N Pseudo-Arc: Finished. Steps=%d, LastParam=%.6f\n', ...
            size(x_branch,2)-1, x_branch(end,end));

    info = build_infoN_local(stop_reason, size(x_branch,2)-1, fail_count_total, arc, wp_arc, n);
end

%% =========================================================
% 加权范数
%% =========================================================
function nrm = weighted_normN_local(v, wp)
    v = v(:);
    n = numel(v) - 1;

    ds = v(1:n);
    dp = v(n+1);

    nrm = sqrt(sum(ds.^2) + (wp^2)*dp^2);
end

%% =========================================================
% 加权内积
%% =========================================================
function val = weighted_dotN_local(a, b, wp)
    a = a(:);
    b = b(:);

    n = numel(a) - 1;

    val = a(1:n).'*b(1:n) + (wp^2)*a(n+1)*b(n+1);
end

%% =========================================================
% 输出 info
%% =========================================================
function info = build_infoN_local(stop_reason, steps, fail_count_total, arc, wp_arc, n)
    info = struct();
    info.stop_reason = stop_reason;
    info.steps = steps;
    info.fail_count_total = fail_count_total;
    info.final_arc = arc;
    info.wp_arc = wp_arc;
    info.state_dim = n;
end