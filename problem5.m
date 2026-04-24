% 初始化方程组参数
A = [1 2 -2 0; 
    1 1  1 0; 
    2 2  1 0; 
    0 0  1 3];
b = [1; 1; 1; 1];
x0 = zeros(4,1); % 初始值零向量
tol = 1e-6;      % 停机条件误差上限
max_iter = 1000; % 最大迭代次数保护

% 计算初始残差范数 ||b - A*x0||，因为x0是零向量，即为 ||b||
norm_b = norm(b - A*x0); 

% 分解矩阵 A = D - L - U
D = diag(diag(A));
L = -tril(A, -1);
U = -triu(A, 1);

%% (1) 雅可比迭代 (Jacobi)
disp('--- 1. 雅可比迭代 (Jacobi) ---');
x_jacobi = x0;
B_J = D \ (L + U);
f_J = D \ b;

for k = 1:max_iter
    x_new = B_J * x_jacobi + f_J;
    if norm(b - A*x_new) / norm_b < tol
        fprintf('Jacobi 迭代收敛！迭代次数: %d\n', k);
        disp('计算结果:'); disp(x_new);
        break;
    end
    x_jacobi = x_new;
end

%% (1) 高斯-赛德尔迭代 (Gauss-Seidel)
disp('--- 2. 高斯-赛德尔迭代 (Gauss-Seidel) ---');
x_gs = x0;
B_GS = (D - L) \ U;
f_GS = (D - L) \ b;

for k = 1:max_iter
    x_new = B_GS * x_gs + f_GS;
    if norm(b - A*x_new) / norm_b < tol
        fprintf('Gauss-Seidel 迭代收敛！迭代次数: %d\n', k);
        disp('计算结果:'); disp(x_new);
        break;
    end
    if norm(x_new) > 1e10 % 发散保护
        fprintf('Gauss-Seidel 迭代发散！已在第 %d 次迭代强行停止。\n', k);
        break;
    end
    x_gs = x_new;
end

%% (2) SOR 迭代
disp('--- 3. SOR 迭代 ---');
% 根据下文的验算，w 必须取极小的欠松弛因子才能收敛，这里选取 w = 0.1
w = 0.1; 
fprintf('选取松弛因子 w = %.1f\n', w);
x_sor = x0;
B_SOR = (D - w*L) \ ((1-w)*D + w*U);
f_SOR = w * ((D - w*L) \ b);

for k = 1:max_iter
    x_new = B_SOR * x_sor + f_SOR;
    if norm(b - A*x_new) / norm_b < tol
        fprintf('SOR 迭代收敛！迭代次数: %d\n', k);
        disp('计算结果:'); disp(x_new);
        break;
    end
    if norm(x_new) > 1e10 % 发散保护
        fprintf('SOR 迭代发散！已在第 %d 次迭代强行停止。\n', k);
        break;
    end
    x_sor = x_new;
end