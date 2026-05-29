%% =========================================================
% Validate_03_AFT_Cubic_Projection.m
%
% 目的：
% 验证三次非线性项 x^3 的 AFT 投影是否正确。
%
% 验证内容：
% 1) x = A cos(t) 的解析投影；
% 2) x = B sin(t) 的解析投影；
% 3) 随机 0/1/3 谐波输入下，与高精度独立投影结果比较。
%
% 谐波系数顺序：
% u = [dc, cos1, sin1, cos3, sin3]^T
% =========================================================

clc; clear; close all;

%% -----------------------------
% Case A: x = A cos(t)
%% -----------------------------
A = 0.7;
u = [0; A; 0; 0; 0];

cubic_aft = cubic_proj_013_local(u);

cubic_theory = [0;
                3*A^3/4;
                0;
                A^3/4;
                0];

err_abs_A = norm(cubic_aft - cubic_theory);
err_rel_A = err_abs_A / max(1e-14, norm(cubic_theory));

fprintf('\n========== Case A: x = A cos(t) ==========\n');
disp('AFT result:');
disp(cubic_aft.');

disp('Theory result:');
disp(cubic_theory.');

fprintf('abs error = %.6e\n', err_abs_A);
fprintf('rel error = %.6e\n', err_rel_A);

%% -----------------------------
% Case B: x = B sin(t)
%% -----------------------------
B = 0.6;
u = [0; 0; B; 0; 0];

cubic_aft = cubic_proj_013_local(u);

cubic_theory = [0;
                0;
                3*B^3/4;
                0;
               -B^3/4];

err_abs_B = norm(cubic_aft - cubic_theory);
err_rel_B = err_abs_B / max(1e-14, norm(cubic_theory));

fprintf('\n========== Case B: x = B sin(t) ==========\n');
disp('AFT result:');
disp(cubic_aft.');

disp('Theory result:');
disp(cubic_theory.');

fprintf('abs error = %.6e\n', err_abs_B);
fprintf('rel error = %.6e\n', err_rel_B);

%% -----------------------------
% Case C: 随机谐波系数 vs 高精度独立投影
%% -----------------------------
rng(2);
u = 0.5 * randn(5,1);

cubic_aft = cubic_proj_013_local(u);
cubic_ref = cubic_proj_013_highres(u, 4096);

err_abs_C = norm(cubic_aft - cubic_ref);
err_rel_C = err_abs_C / max(1e-14, norm(cubic_ref));

fprintf('\n========== Case C: random u vs high-res projection ==========\n');
disp('u:');
disp(u.');

disp('AFT result:');
disp(cubic_aft.');

disp('High-res reference:');
disp(cubic_ref.');

fprintf('abs error = %.6e\n', err_abs_C);
fprintf('rel error = %.6e\n', err_rel_C);

%% -----------------------------
% 判定
%% -----------------------------
fprintf('\n========== 判定 ==========\n');

if err_rel_A < 1e-12
    fprintf('Case A 通过：cos 基波三次项投影正确。\n');
else
    fprintf('Case A 未通过：请检查 cos1/cos3 系数或投影归一化。\n');
end

if err_rel_B < 1e-12
    fprintf('Case B 通过：sin 基波三次项投影正确。\n');
else
    fprintf('Case B 未通过：请检查 sin1/sin3 符号。\n');
end

if err_rel_C < 1e-10
    fprintf('Case C 通过：一般谐波输入下 AFT 投影与高精度投影一致。\n');
else
    fprintf('Case C 未通过：请检查 T_mat/T_inv 或采样点数。\n');
end

fprintf('\n========== AFT 三次项投影验证完成 ==========\n');

%% =========================================================
% 本脚本局部函数
%% =========================================================

function cubic = cubic_proj_013_local(u)
    N = 64;
    t = (0:N-1)'*(2*pi/N);

    dc = ones(N,1);
    c1 = cos(t);
    s1 = sin(t);
    c3 = cos(3*t);
    s3 = sin(3*t);

    T = [dc, c1, s1, c3, s3];

    Inv = [dc, 2*c1, 2*s1, 2*c3, 2*s3]';
    Tinv = (1/N) * Inv;
    Tinv(1,:) = (1/N) * dc';

    x_time = T * u;
    cubic = Tinv * (x_time.^3);
end

function cubic = cubic_proj_013_highres(u, N)
    t = (0:N-1)'*(2*pi/N);

    dc = ones(N,1);
    c1 = cos(t);
    s1 = sin(t);
    c3 = cos(3*t);
    s3 = sin(3*t);

    T = [dc, c1, s1, c3, s3];

    x_time = T * u;
    x3 = x_time.^3;

    cubic = zeros(5,1);
    cubic(1) = mean(x3);
    cubic(2) = 2*mean(x3 .* c1);
    cubic(3) = 2*mean(x3 .* s1);
    cubic(4) = 2*mean(x3 .* c3);
    cubic(5) = 2*mean(x3 .* s3);
end