% -*- coding: gbk -*-
% File          : RIV.m
% Creation Date : 2015-05-07
% Description   : 辅助变量最小二乘法
% 

for k = nMax + 1 : L + nMax
    for i = 1:na
        h(i, k) = -z(k - i);
        h1(i, k) = -x(k - i);
    end
    
    for i = 1:nb
        h(na + i, k) = u(k - i);
        h1(na + i, k) = u(k - i);
    end
    
    s(k) = h(:, k)' * P(:, :, k-1) * h(:, k) + 1.0;
    Inn(k) = z(k) - h(:, k)' * Theta(:, k-1);
    K(:, k) = P(:, :, k-1) * h1(:, k) / s(k);
    P(:, :, k) = P(:, :, k-1) - K(:, k) * h(:, k)' * P(:, :, k-1);
    Theta(:, k) = Theta(:, k-1) + K(:, k) * Inn(k);
    J(k) = J(k-1) + Inn(k)^2 / s(k);
    
    if k > dTime
        ThetaIV(:, k) = (1 - alpha) * ThetaIV(:, k-1) + alpha * Theta(:, k - dTime);
    else
        ThetaIV(:, k) = Theta(:, k);
    end
    
    x(k) = h1(:, k)' * ThetaIV(:, k);
end
