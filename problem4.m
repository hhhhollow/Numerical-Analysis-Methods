load('data2.mat');
A_100 = full(A(1:100, 1:100)); 
B_100 = B(1:100, 1);
n = 100;


%% 平方根法
L = zeros(n, n);

% 计算 L 矩阵
for j = 1:n
    % （3.20）l_jj
    sum_sq = 0;
    for k = 1:j-1
        sum_sq = sum_sq + L(j, k)^2;
    end
    L(j, j) = sqrt(A_100(j, j) - sum_sq);
    
    % （3.21）l_ij
    for i = j+1:n
        sum_cross = 0;
        for k = 1:j-1
            sum_cross = sum_cross + L(i, k) * L(j, k);
        end
        L(i, j) = (A_100(i, j) - sum_cross) / L(j, j);
    end
end



% 解 Ly = b
y1 = zeros(n, 1);
for i = 1:n
    sum_Ly = 0;
    for k = 1:i-1
        sum_Ly = sum_Ly + L(i, k) * y1(k);
    end
    y1(i) = (B_100(i) - sum_Ly) / L(i, i);
end

% 解 L^T x = y
x1 = zeros(n, 1);
for i = n:-1:1
    sum_Ltx = 0;
    for k = i+1:n
        sum_Ltx = sum_Ltx + L(k, i) * x1(k);
    end
    x1(i) = (y1(i) - sum_Ltx) / L(i, i);
end


%% 改进的平方根法
L2 = eye(n); % 主对角线全为 1
D = zeros(n, 1);

% 计算 L 和 D
for j = 1:n
    % 先计算辅助项 t_jk = l_jk * d_k
    t = zeros(j-1, 1);
    for k = 1:j-1
        t(k) = L2(j, k) * D(k);
    end
    
    % 计算 d_j
    sum_lt_d = 0;
    for k = 1:j-1
        sum_lt_d = sum_lt_d + L2(j, k) * t(k);
    end
    D(j) = A_100(j, j) - sum_lt_d;
    
    % 计算 l_ij
    for i = j+1:n
        sum_lt_l = 0;
        for k = 1:j-1
            sum_lt_l = sum_lt_l + L2(i, k) * t(k);
        end
        L2(i, j) = (A_100(i, j) - sum_lt_l) / D(j);
    end
end

% 解 Ly = b
y2 = zeros(n, 1);
for i = 1:n
    sum_Ly2 = 0;
    for k = 1:i-1
        sum_Ly2 = sum_Ly2 + L2(i, k) * y2(k);
    end
    y2(i) = B_100(i) - sum_Ly2;
end

% 解 DL^T x = y
x2 = zeros(n, 1);
for i = n:-1:1
    sum_Ltx2 = 0;
    for k = i+1:n
        sum_Ltx2 = sum_Ltx2 + L2(k, i) * x2(k);
    end
    x2(i) = y2(i) / D(i) - sum_Ltx2;
end


%% ================= 严谨验算环节 =================
% res1 = max(abs(A_100 * x1 - B_100));
% res2 = max(abs(A_100 * x2 - B_100));
% fprintf('平方根法无穷大范数残差: %e\n', res1);
% fprintf('改进法无穷大范数残差: %e\n', res2);