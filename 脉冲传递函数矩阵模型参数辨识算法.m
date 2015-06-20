% -*- coding: gbk -*-
% File          : 脉冲传递函数矩阵模型参数辨识算法.m
% Creation Date : 2015-06-20
% Description   : 脉冲传递函数矩阵模型参数辨识算法, chapter 12
% 

for k = 1+n : L+n
    for i = 1:m
        % 初始化
        for p = 1:N
            h(p, k) = 0.0;
        end
        
        % 构造数据向量
        for p = 1:n
            h(p, k) = -z(i, k-p);
        end
        for j = 1:r
            for p = 1:n
                h(n + (i-1) * r * n + (j-1) * n + p, k) = u(j, k-p);
            end
        end
        
        % 辨识算法
        if i == 1
            P(:, :, k) = P(:, :, k-1);
            Theta(:, k) = Theta(:, k-1);
        end
        s(k) = h(:, k)' * P(:, :, k) * h(:, k) + 1.0;
        Inn(k) = z(i, k) - h(:, k)' * Theta(:, k);
        K(:, k) = P(:, :, k) * h(:, k) / s(k);
        P(:, :, k) = P(:, :, k) - K(:, k) * K(:, k)' * s(k);
        Theta(:, k) = Theta(:, k) + K(:, k) * Inn(k);
        
        % 第i子系统损失函数和各噪声方差及各子系统噪声方差和
        J(i, k) = J(i, k-1) + Inn(k)^2 / s(k);
        sigma(i, k) = J(i, k) / (1 + ThetaC(:, k-1)' * D * Theta(:, k));
        if i == m
            Sigma = 0.0;
            for j = 1:m
                Sigma = Sigma + sigma(j, k);
            end
            ThetaC(:, k) = Theta(:, k) + Sigma * P(:, :, k) * D * Theta(:, k-1);
        end
    end
end
