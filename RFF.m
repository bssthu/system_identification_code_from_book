% -*- coding: gbk -*-
% File          : RFF.m
% Creation Date : 2015-06-12
% Description   : “≈Õ¸“Ú◊”∑®, chapter 5
% 

for k = nMax + 1 : L + nMax
    for i = 1:na
        h(i, k) = -z(k - i);
    end
    
    for i = 1:nb
        h(na + i, k) = u(k - i);
    end
    
    s(k) = h(:, k)' * P(:, :, k-1) * h(:, k) + mu;
    Inn(k) = z(k) - h(:, k)' * Theta(:, k-1);
    K(:, k) = P(:, :, k-1) * h(:, k) / s(k);
    P(:, :, k) = (P(:, :, k-1) - K(:, k) * K(:, k)' * s(k)) / mu;
    Theta(:, k) = Theta(:, k-1) + K(:, k) * Inn(k);
    
    J(k) = J(k-1) + Inn(k)^2 / s(k);
end
