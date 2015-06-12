% -*- coding: gbk -*-
% File          : RELS.m
% Creation Date : 2015-05-07
% Description   : 增广最小二乘法, chapter 6
% 

for k = nMax + 1 : L + nMax
    % 构造数据向量
    for i = 1:na
        h(i, k) = -z(k - i);
    end
    for i = 1:nb
        h(na + i, k) = u(k - i);
    end
    for i = 1:nd
        h(na + nb + i, k) = v1(k - i);
    end
    
    % 辨识算法
    s(k) = h(:, k)' * P(:, :, k-1) * h(:, k) + 1.0;
    Inn(k) = z(k) - h(:, k)' * Theta(:, k-1);
    K(:, k) = P(:, :, k-1) * h(:, k) / s(k);
    P(:, :, k) = P(:, :, k-1) - K(:, k) * K(:, k)' * s(k);
    Theta(:, k) = Theta(:, k-1) + K(:, k) * Inn(k);
    
    % 损失函数
    J(k) = J(k-1) + Inn(k)^2 / s(k);
    
    v1(k) = z(k) - h(:, k)' * Theta(:, k);
end
