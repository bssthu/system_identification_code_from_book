% -*- coding: gbk -*-
% File          : RFF.m
% Creation Date : 2015-06-12
% Description   : 遗忘因子法, chapter 5
% 

for k = nMax + 1 : L + nMax
    % 构造数据向量
    for i = 1:na
        h(i, k) = -z(k - i);
    end
    
    for i = 1:nb
        h(na + i, k) = u(k - i);
    end
    
    % 辨识算法
    s(k) = h(:, k)' * P(:, :, k-1) * h(:, k) + mu;
    Inn(k) = z(k) - h(:, k)' * Theta(:, k-1);
    K(:, k) = P(:, :, k-1) * h(:, k) / s(k);
    P(:, :, k) = (P(:, :, k-1) - K(:, k) * K(:, k)' * s(k)) / mu;
    Theta(:, k) = Theta(:, k-1) + K(:, k) * Inn(k);
    
    % 损失函数
    J(k) = J(k-1) + Inn(k)^2 / s(k);
end
