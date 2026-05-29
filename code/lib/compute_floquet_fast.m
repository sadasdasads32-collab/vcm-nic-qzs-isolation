function maxMu = compute_floquet_fast(x_coeff, sysP, Omega, Nt_floquet)
% 快速 Floquet 乘子计算：基于 HB 轨道重构 A(t)，RK4 积分单周期矩阵
    x1c = x_coeff(1:5);
    x2c = x_coeff(6:10);
    qc  = x_coeff(11:15);

    Tp = 2*pi/Omega;
    dt = Tp / Nt_floquet;

    Phi = eye(6);
    t = 0;
    for i = 1:Nt_floquet
        A = build_A_matrix(t, sysP, Omega, x1c, x2c, qc);
        k1 = A * Phi;
        A_mid = build_A_matrix(t+0.5*dt, sysP, Omega, x1c, x2c, qc);
        k2 = A_mid * (Phi + 0.5*dt*k1);
        k3 = A_mid * (Phi + 0.5*dt*k2);
        A_end = build_A_matrix(t+dt, sysP, Omega, x1c, x2c, qc);
        k4 = A_end * (Phi + dt*k3);
        Phi = Phi + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
        t = t + dt;
    end
    maxMu = max(abs(eig(Phi)));
end

function A = build_A_matrix(t, sysP, Omega, x1c, x2c, qc)
    w = Omega;
    ct=cos(w*t); st=sin(w*t); c3t=cos(3*w*t); s3t=sin(3*w*t);

    x1 = x1c(1)+x1c(2)*ct+x1c(3)*st+x1c(4)*c3t+x1c(5)*s3t;
    x2 = x2c(1)+x2c(2)*ct+x2c(3)*st+x2c(4)*c3t+x2c(5)*s3t;

    be1=sysP(1); be2=sysP(2); mu=sysP(3);
    al1=sysP(4); ga1=sysP(5); ze=sysP(6);
    lam=sysP(7); kap_e=sysP(8); kap_c=sysP(9); sigma=sysP(10); ga2=sysP(11);

    theta = sqrt(max(lam,0));
    dx = x1-x2;

    df12 = (be1+al1) + 3*ga1*dx^2;
    df2g_x = be2 + 3*ga2*x2^2;
    df2g_v = 2*mu*ze;

    A = zeros(6);
    A(1,2) = 1;
    A(2,1) = -df12;       A(2,3) = df12;                        A(2,6) = theta;
    A(3,4) = 1;
    A(4,1) = df12/mu;     A(4,3) = (-df12 - df2g_x)/mu;         A(4,4) = -df2g_v/mu;  A(4,6) = -theta/mu;
    A(5,6) = 1;

    if abs(kap_e) > 1e-14
        A(6,2) = -theta/kap_e;   A(6,4) = theta/kap_e;
        A(6,5) = -kap_c/kap_e;   A(6,6) = -sigma/kap_e;
    else
        s = max(abs(sigma), 1e-12);
        A(6,2) = -50*theta/s;    A(6,4) = 50*theta/s;
        A(6,5) = -50*kap_c/s;    A(6,6) = -50;
    end
end
