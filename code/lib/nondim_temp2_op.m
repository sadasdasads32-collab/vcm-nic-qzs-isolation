function [z, Jac] = nondim_temp2_op(y, sysP)
%% nondim_temp2_op
% 消元后的 10 维复算子模型残差
%
% 输入:
%   y    : 11x1 = [x1(5); x2(5); Omega_or_Fw]
%   sysP : 11x1，与 15 维完整模型相同
%
% 输出:
%   z    : 10x1 机械残差 [R1; R2]
%   Jac  : 10x10 解析 Jacobian
%
% 与当前 15 维完整模型保持一致的电路方程:
%
%   kap_e*Q'' + sigma*Q' + kap_c*Q - theta*(x1' - x2') = 0
%
% 系数域写作:
%
%   A*q - theta*D*x12 = 0
%
% 因此:
%
%   q = theta*A^{-1}*D*x12
%
% 电磁反馈力:
%
%   force_em = theta*q'
%            = theta*D*q
%            = theta^2*D*A^{-1}*D*x12
%
% 注意:
%   这里的符号必须和 nondim_temp2.m 中的 R3 保持一致。
%   如果 15 维模型中 R3 = A*q - theta*D*x12，
%   那么这里必须使用 q = +theta*A^{-1}*D*x12。

    global Fw FixedOmega

    %% =========================
    % 1) 输入检查与模式切换
    %% =========================
    if numel(y) ~= 11
        error('nondim_temp2_op expects y to be 11x1: [x1(5); x2(5); Omega_or_Fw].');
    end

    if numel(sysP) ~= 11
        error('nondim_temp2_op expects sysP to be 11x1.');
    end

    state = y(1:10);

    if isempty(FixedOmega)
        W = y(11);          % 扫频模式：y(11)=Omega
        current_Fw = Fw;    % 激励幅值固定
    else
        W = FixedOmega;     % 扫力模式：Omega固定
        current_Fw = y(11); % y(11)=Fw
    end

    %% =========================
    % 2) 参数映射
    %% =========================
    be1 = sysP(1);
    be2 = sysP(2);
    mu_mass = sysP(3);

    al1 = sysP(4);
    ga1 = sysP(5);
    ze1 = sysP(6);

    lam_phys = sysP(7);
    kap_e = sysP(8);
    kap_c = sysP(9);
    sigma = sysP(10);
    ga2 = sysP(11);

    theta = sqrt(max(0, lam_phys));

    %% =========================
    % 2.5) 断耦合基准阻尼
    %% =========================
    % 与 nondim_temp2.m 保持一致:
    % lam > 0 时不启用额外层间阻尼；
    % lam = 0 时启用断耦合基准层间阻尼。
    lam_eps = 1e-12;
    zeta12 = 0.0;

    if abs(lam_phys) < lam_eps
        zeta12 = 0.05;
        % 如果你的 15维 nondim_temp2.m 中这里使用 0.020412
        % 或自动换算公式，则这里也必须同步修改。
    end

    %% =========================
    % 3) 状态拆分
    %% =========================
    x1 = state(1:5);
    x2 = state(6:10);

    x12 = x1 - x2;

    cubic12 = cubic_proj_013(x12);
    cubic2  = cubic_proj_013(x2);

    %% =========================
    % 4) 谐波导数算子
    %% =========================
    W2 = W^2;

    % 二阶导数算子
    Mat_Inertia = diag([0; -W2; -W2; -9*W2; -9*W2]);

    % 一阶导数算子
    Mat_Deriv = zeros(5);
    Mat_Deriv(2,3) = W;
    Mat_Deriv(3,2) = -W;
    Mat_Deriv(4,5) = 3*W;
    Mat_Deriv(5,4) = -3*W;

    I5 = eye(5);

    %% =========================
    % 5) 电路自由度严格消元
    %% =========================
    % 当前 15维完整模型的电路方程为:
    %
    %   kap_e*q'' + sigma*q' + kap_c*q - theta*x12' = 0
    %
    % 系数域:
    %
    %   A*q - theta*D*x12 = 0
    %
    % 因此:
    %
    %   q = theta*A^{-1}*D*x12

    A = kap_e*Mat_Inertia + sigma*Mat_Deriv + kap_c*I5;

    Dx12 = Mat_Deriv * x12;

    q = theta * (A \ Dx12);

    q_dot = Mat_Deriv * q;

    % 机械端电磁反馈力
    force_em = theta * q_dot;

    %% =========================
    % 6) 阻尼项
    %% =========================
    x12_dot = Mat_Deriv * x12;
    damp12  = (2*zeta12) * x12_dot;

    x2_dot = Mat_Deriv * x2;
    damp2  = (2*mu_mass*ze1) * x2_dot;

    %% =========================
    % 7) 机械残差 R1, R2
    %% =========================
    % R1:
    % x1'' + (be1+al1)*x12 + ga1*x12^3
    %      + theta*q' + damp12 = Fw*cos

    R1 = (Mat_Inertia*x1) ...
         + (be1+al1)*x12 ...
         + ga1*cubic12 ...
         + force_em ...
         + damp12;

    R1(2) = R1(2) - current_Fw;

    % R2:
    % mu*x2'' + damp2 + be2*x2 + ga2*x2^3
    % - [(be1+al1)*x12 + ga1*x12^3]
    % - theta*q' - damp12 = 0

    Force_from_upper = (be1+al1)*x12 + ga1*cubic12;

    R2 = mu_mass*(Mat_Inertia*x2) ...
         + damp2 ...
         + be2*x2 ...
         + ga2*cubic2 ...
         - Force_from_upper ...
         - force_em ...
         - damp12;

    z = [R1; R2];

    %% =========================
    % 8) 解析 Jacobian
    %% =========================
    if nargout > 1

        J_cubic_x12 = AFT_GetJac(x12);
        J_cubic_x2  = AFT_GetJac(x2);

        % force_em = theta^2 * D * A^{-1} * D * (x1 - x2)
        %
        % 定义:
        %   Op = theta^2 * D * A^{-1} * D
        %
        % 则:
        %   d(force_em)/dx1 =  Op
        %   d(force_em)/dx2 = -Op

        Op = theta^2 * (Mat_Deriv * (A \ Mat_Deriv));

        Jf_x1 = Op;
        Jf_x2 = -Op;

        %% ---- R1 Jacobian ----
        J11 = Mat_Inertia ...
              + (be1+al1)*I5 ...
              + ga1*J_cubic_x12 ...
              + Jf_x1;

        J12 = -(be1+al1)*I5 ...
              - ga1*J_cubic_x12 ...
              + Jf_x2;

        if zeta12 ~= 0
            J11 = J11 + (2*zeta12)*Mat_Deriv;
            J12 = J12 - (2*zeta12)*Mat_Deriv;
        end

        %% ---- R2 Jacobian ----
        % R2 中含 -force_em
        % 所以对 x1 的导数是 -Jf_x1
        % 对 x2 的导数是 -Jf_x2

        J21 = -(be1+al1)*I5 ...
              - ga1*J_cubic_x12 ...
              - Jf_x1;

        J22 = mu_mass*Mat_Inertia ...
              + (2*mu_mass*ze1)*Mat_Deriv ...
              + be2*I5 ...
              + ga2*J_cubic_x2 ...
              + (be1+al1)*I5 ...
              + ga1*J_cubic_x12 ...
              - Jf_x2;

        if zeta12 ~= 0
            J21 = J21 - (2*zeta12)*Mat_Deriv;
            J22 = J22 + (2*zeta12)*Mat_Deriv;
        end

        Jac = [J11, J12;
               J21, J22];
    end
end

%% =========================================================
% AFT: 三次项投影
%% =========================================================
function cubic = cubic_proj_013(u)

    [~, T_mat, T_inv] = get_AFT_matrices();

    cubic = T_inv * ((T_mat * u).^3);
end

%% =========================================================
% AFT: 三次项 Jacobian
%% =========================================================
function J_aft = AFT_GetJac(u)

    [~, T_mat, T_inv] = get_AFT_matrices();

    u_time = T_mat * u;

    df_du = 3 * u_time.^2;

    % 等价写法:
    % J_aft = T_inv * diag(df_du) * T_mat;
    % 为了和你的 nondim_temp2.m 写法风格一致，也可以写成下面这样。
    J_aft = T_inv * (df_du .* T_mat);
end

%% =========================================================
% AFT 矩阵
%% =========================================================
function [N, T_mat, T_inv] = get_AFT_matrices()

    persistent pN pT pTinv

    if isempty(pN)

        pN = 64;

        t = (0:pN-1)'*(2*pi/pN);

        dc = ones(pN,1);

        c1 = cos(t);
        s1 = sin(t);

        c3 = cos(3*t);
        s3 = sin(3*t);

        pT = [dc, c1, s1, c3, s3];

        Inv = [dc, 2*c1, 2*s1, 2*c3, 2*s3]';

        pTinv = (1/pN) * Inv;
        pTinv(1,:) = (1/pN) * dc';
    end

    N = pN;
    T_mat = pT;
    T_inv = pTinv;
end