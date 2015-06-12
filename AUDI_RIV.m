% -*- coding: gbk -*-
% File          : AUDI_RIV.m
% Creation Date : 2015-06-12
% Description   : AUDI 辅助变量最小二乘法, chapter 11
% 

for k = 1+n : L+n
    % 构造增广数据向量
    for i = 0:n-1
        Phi(2*i + 1, k) = -z(k - n + i);
        Phi1(2*i + 1, k) = -x(k - n + i);
        h1(2*i + 1, k) = -x(k - n + i);
        Phi(2*i + 2, k) = u(k - n + i);
        Phi1(2*i + 2, k) = -x(k - n + i);
        h1(2*i + 2, k) = u(k - n + i);
    end
    Phi(2*n + 1, k) = -z(k);
    
    % 计算辅助变量
    if k >= k0
        for i = 1:N-1
            x(k) = x(k) + (alpha * U(i, N, k-1) + ...
                    (1 - alpha) * U(i, N, k-1)) * h1(i, k);
        end
        if abs(x(k)) >= 1.1 * abs(z(k))
            x(k) = z(k);
        end
    else
        x(k) = z(k);
    end
    
    % AUDI辨识算法
    Phi1(2*n + 1, k) = -x(k);
    f(:, k) = U(:, :, k-1)' * Phi(:, k);
    g(:, k) = D(:, :, k-1) * f(:, k);
    f1(:, k) = V(:, :, k-1)' * Phi1(:, k);
    g1(:, k) = D(:, :, k-1) * f1(:, k);
    Beta(1) = 1.0;
    
    for j = 1:N
        Beta(j+1) = Beta(j) + f(j, k) * g1(j, k);
        D(j, j, k) = D(j, j, k-1) * Beta(j) / Beta(j+1);
        E(j) = -f(j, k) / Beta(j);
        E1(j) = -f1(j, k) / Beta(j);
        K(j) = g1(j, k);
        K1(j) = g(j, k);
        
        for i = 1:j-1
            U(i, j, k) = U(i, j, k-1) + K(i) * E(j);
            K(i) = K(i) + U(i, j, k-1) * g1(j, k);
            V(i, j, k) = V(i, j, k-1) + K1(i) * E1(j);
            K1(i) = K1(i) + V(i, j, k-1) * g(j, k);
        end
        
        U(j, j, k) = 1.0;
        V(j, j, k) = 1.0;
    end
    
    % 各阶模型损失函数
    for r = 1:n
        J(r, k) = J(r, k-1) + (f(2*r + 1, k)^2) / Beta(2*r + 1);
    end
end

% 各阶模型噪声标准差估计
for r = 1:n
    LambdaJ(r) = sqrt(J(r, L+n) / L);
end
