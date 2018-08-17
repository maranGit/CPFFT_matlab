function C2 = ddot42_2(A4, B2, N)

vecflag = 0;
if isvector(B2)
    vecflag = 1;
    B2 = reshape(B2, N, 9);
end

C2 = zeros(N,9);
for ii = 1:N
    C2(ii,:) = B2(ii,:) * reshape(A4(ii,:),9,9);
end

if vecflag
    C2 = reshape(C2, N*9, 1);
end

end