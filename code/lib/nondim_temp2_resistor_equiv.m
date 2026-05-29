function [z, Jac] = nondim_temp2_resistor_equiv(y, sysP)
%% nondim_temp2_resistor_equiv
% 纯电阻电路的理论等效电磁阻尼模型
%
% 目的：
%   用来验证 nondim_temp2_1 中“纯电阻电路支路”是否正确。
%
% 理论关系：
%   纯电阻电路方程：
%       sigma*q_dot - theta*(x1_dot - x2_dot) = 0
%
%   因此：
%       q_dot = theta/sigma * (x1_dot - x2_dot)
%
%   机械方程中的电磁力：
%       theta*q_dot = theta^2/sigma * (x1_dot - x2_dot)
%
%   因为：
%       theta^2 = lambda
%
%   所以纯电阻电路等效为层间速度阻尼：
%       c_em = lambda/sigma
%
%   若写成：
%       2*zeta_em*(x1_dot - x2_dot)
%
%   则：
%       zeta_em = lambda/(2*sigma)
%
% 输入:
%   y    : 16x1 = [x1(5); x2(5); q(5); Omega_or_Fw]
%   sysP : 11x1 参数向量
%
% 参数排列:
%   sysP = [be1, be2, mu, al1, ga1, ze1, lam, kap_e, kap_c, sigma, ga2]
%
% 输出:
%   z    : 15x1 残差
%   Jac  : 15x15 雅可比
%
% 说明：
%   本函数没有真正使用电路变量 q。
%   为了保持和原程序相同的 15 维未知量结构，令：
%       R3 = q
%   即强制 q = 0。
%
%   这样可以继续直接调用：
%       newton('nondim_temp2_resistor_equiv', ...)
%       branch_follow2('nondim_temp2_resistor_equiv', ...)

    global Fw FixedOmega

%% =========================
% 1) 输入检查 & 模式切换
% =========================
    if numel(y) ~= 16
        error('nondim_temp2_resistor_equiv expects y to be 16x1: [x1(5); x2(5); q(5); Omega_or_Fw].');
    end

    if numel(sysP) ~= 11
        error('nondim_temp2_resistor_equiv expects sysP to be 11x1.');
    end

    state = y(1:15);

    if isempty(FixedOmega)
        W = y(16);          % 扫频：y(16)=Omega
        current_Fw = Fw;    % 激励幅值固定
    else
        W = FixedOmega;     % 扫力：Omega固定
        current_Fw = y(16); % y(16)=Fw
    end

%% =========================
% 2) 参数映射
% =========================
    be1     = sysP(1);
    be2     = sysP(2);
    mu_mass = sysP(3);

    al1 = sysP(4);
    ga1 = sysP(5);
    ze1 = sysP(6);

    lam_phys = sysP(7);
    sigma    = sysP(10);
    ga2      = sysP(11);

    if abs(sigma) < 1e-12
        error('sigma is zero. Equivalent resistor damping requires sigma > 0.');
    end

    % 纯电阻等效层间阻尼：
    % c_em = lambda/sigma = 2*zeta_em
    c_em = lam_phys/sigma;

%% =========================
% 3) 状态拆分 & 基础算子
% =========================
    x1 = state(1:5);
    x2 = state(6:10);
    q  = state(11:15);

    x12 = x1 - x2;

    cubic12 = cubic_proj_013(x12);
    cubic2  = cubic_proj_013(x2);

    W2 = W^2;

    % 一阶导数算子：对 [a0,a1,b1,a3,b3]
    Mat_Deriv = zeros(5);
    Mat_Deriv(2,3) = W;    Mat_Deriv(3,2) = -W;
    Mat_Deriv(4,5) = 3*W;  Mat_Deriv(5,4) = -3*W;

    % 二阶导数算子
    Mat_Inertia = diag([0; -W2; -W2; -9*W2; -9*W2]);

    I5 = eye(5);

%% =========================
% 4) 等效阻尼项
% =========================
    x12_dot = Mat_Deriv*x12;
    x2_dot  = Mat_Deriv*x2;

    % 等效电磁阻尼力
    damp_em = c_em*x12_dot;

    % 下层对地阻尼
    damp2 = (2*mu_mass*ze1)*x2_dot;

    % 层间非线性恢复力
    Force_from_upper = (be1+al1)*x12 + ga1*cubic12;

%% =========================
% 5) 残差 R1, R2, R3
% =========================

    % 上层方程：
    % x1'' + Force_from_upper + c_em*(x1' - x2') = Fw*cos
    R1 = Mat_Inertia*x1 ...
         + Force_from_upper ...
         + damp_em;

    R1(2) = R1(2) - current_Fw;

    % 下层方程：
    % mu*x2'' + 2*mu*ze1*x2' + be2*x2 + ga2*x2^3
    % - Force_from_upper - c_em*(x1' - x2') = 0
    R2 = mu_mass*(Mat_Inertia*x2) ...
         + damp2 ...
         + be2*x2 + ga2*cubic2 ...
         - Force_from_upper ...
         - damp_em;

    % q 在等效模型中不参与动力学。
    % 为保持 15 维结构，强制 q = 0。
    R3 = q;

    z = [R1; R2; R3];

%% =========================
% 6) 解析雅可比
% =========================
    if nargout > 1
        J_cubic_x12 = AFT_GetJac(x12);
        J_cubic_x2  = AFT_GetJac(x2);

        % R1
        J11 = Mat_Inertia ...
              + (be1+al1)*I5 ...
              + ga1*J_cubic_x12 ...
              + c_em*Mat_Deriv;

        J12 = -(be1+al1)*I5 ...
              - ga1*J_cubic_x12 ...
              - c_em*Mat_Deriv;

        J13 = zeros(5);

        % R2
        J21 = -(be1+al1)*I5 ...
              - ga1*J_cubic_x12 ...
              - c_em*Mat_Deriv;

        J22 = mu_mass*Mat_Inertia ...
              + (2*mu_mass*ze1)*Mat_Deriv ...
              + be2*I5 ...
              + ga2*J_cubic_x2 ...
              + (be1+al1)*I5 ...
              + ga1*J_cubic_x12 ...
              + c_em*Mat_Deriv;

        J23 = zeros(5);

        % R3 = q
        J31 = zeros(5);
        J32 = zeros(5);
        J33 = eye(5);

        Jac = [J11, J12, J13;
               J21, J22, J23;
               J31, J32, J33];
    end
end

%% ============================================================
% AFT 辅助函数
% ============================================================
function cubic = cubic_proj_013(u)
    [~, T_mat, T_inv] = get_AFT_matrices();
    cubic = T_inv * ((T_mat * u).^3);
end

function J_aft = AFT_GetJac(u)
    [~, T_mat, T_inv] = get_AFT_matrices();
    u_time = T_mat * u;
    df_du = 3 * u_time.^2;
    J_aft = T_inv * (df_du .* T_mat);
end

function [N, T_mat, T_inv] = get_AFT_matrices()
    persistent pN pT pTinv

    if isempty(pN)
        pN = 64;
        t = (0:pN-1)'*(2*pi/pN);

        c1 = cos(t);
        s1 = sin(t);
        c3 = cos(3*t);
        s3 = sin(3*t);
        dc = ones(pN,1);

        pT = [dc, c1, s1, c3, s3];

        Inv = [dc, 2*c1, 2*s1, 2*c3, 2*s3]';
        pTinv = (1/pN) * Inv;
        pTinv(1,:) = (1/pN) * dc';
    end

    N = pN;
    T_mat = pT;
    T_inv = pTinv;
end
