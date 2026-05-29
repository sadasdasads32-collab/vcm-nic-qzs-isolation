function TF = compute_TF_fast(x_coeff, sysP, Omega, Fw_val)
% 快速力传递率计算（仅基波，与 optimization.m 一致）
    Nt = 512;
    x2c = x_coeff(6:10);
    be2 = sysP(2);
    mu  = sysP(3);
    ze  = sysP(6);
    ga2 = sysP(11);

    Tp = 2*pi/Omega;
    t = linspace(0, Tp, Nt+1); t(end) = [];

    ct = cos(Omega*t);  st = sin(Omega*t);
    c3t = cos(3*Omega*t); s3t = sin(3*Omega*t);

    x2  = x2c(1) + x2c(2)*ct + x2c(3)*st + x2c(4)*c3t + x2c(5)*s3t;
    v2  = -Omega*x2c(2)*st + Omega*x2c(3)*ct ...
          -3*Omega*x2c(4)*s3t + 3*Omega*x2c(5)*c3t;

    Ftr = be2*x2 + ga2*(x2.^3) + 2*mu*ze*v2;
    Y = fft(Ftr)/Nt;
    A1 = 2*abs(Y(2));
    TF = A1 / max(Fw_val, 1e-12);
end
