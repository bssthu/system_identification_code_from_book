% -*- coding: gbk -*-
% File          : RGLS.m
% Creation Date : 2015-06-12
% Description   : 广义最小二乘法, chapter 6
% 

for k = nMax + 1 : L + nMax
    % 输入、输出滤波值
    zf(k) = z(k);
    uf(k) = u(k);
    for i = 1:nc
        zf(k) = zf(k) + Thetae(i, k - 1) * z(k-i);
        uf(k) = uf(k) + Thetae(i, k - 1) * u(k-i);
    end
    
    %构造数据向量和滤波数据向量
    for i = 1:na
        hf(i, k) = -zf(k - i);
        h(i, k) = -z(k - i);
    end
    
    for i = 1:nb
        hf(na + i, k) = uf(k - i);
        h(na + i, k) = u(k - i);
    end
    
    for i = 1:nc
        he(i, k) = -e1(k - i);
    end
    
    % 辨识算法
    s(k) = hf(:, k)' * P(:, :, k-1) * hf(:, k) + 1.0;
    Inn(k) = zf(k) - hf(:, k)' * Theta(:, k-1);
    K(:, k) = P(:, :, k-1) * hf(:, k) / s(k);
    P(:, :, k) = P(:, :, k-1) - K(:, k) * K(:, k)' * s(k);
    Theta(:, k) = Theta(:, k-1) + K(:, k) * Inn(k);
    se(k) = he(:, k)' * Pe(:, :, k-1) * he(:, k) + 1.0;
    e1(k) = z(k) - h(:, k)' * Theta(:, k);
    Inne(k) = e1(k) - he(:, k)' * Thetae(:, k - 1);
    Ke(:, k) = Pe(:, :, k - 1) * he(:, k) / se(k);
    Pe(:, :, k) = Pe(:, :, k - 1) - Ke(:, k) * Ke(:, k)' * se(k);
    Thetae(:, k) = Thetae(:, k - 1) + Ke(:, k) * Inne(k);
    
    % 损失函数
    J(k) = J(k-1) + Inn(k)^2 / s(k);
    Je(k) = Je(k-1) + Inne(k)^2 / se(k);
end
