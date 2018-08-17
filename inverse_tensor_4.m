function K4_inv = inverse_tensor_4(K4, n)
K4_inv = zeros(n,81);
for ii = 1:n
    K4_99 = reshape(K4(ii,:),9,9);
    K4_inv_99 = inv(K4_99);
    K4_inv(ii,:) = reshape(K4_inv_99,1,81);
end
end