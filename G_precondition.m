function G_K_dFvector = G_precondition(A2vector, Ghat4, N, ndim, K4, K_pre, K_pre_inv)
%
% problem size
ndim2 = ndim ^ 2;
N3 = N ^ 3;

% =========================================================================
% ========================== precondition =================================
% =========================================================================
% 
% A2vector = ddot42_2( K_pre_inv, A2vector, N3 );
A2 = reshape( A2vector, N3, ndim2 );
A2 = ddot42_2( K_pre_inv, A2, N3 );
% for temp = 1:ndim2
%     A2(:,temp) = A2(:,temp) .* cg_precondition;
% end
% A2 = ddot42_2(cg_precondition_inv,A2, N3);
% =========================================================================
% =========================================================================
% =========================================================================
B2 = ddot42_2( K4, A2, N3 );

% (inverse) Fourier transform (for each tensor component in each direction)
fftfem = @(x) reshape( fftshift( fftn( ifftshift( reshape( x, N, N, N ) ) ) ), N3, 1 );
ifftfem = @(x) reshape( fftshift( ifftn( ifftshift( reshape( x, N, N, N ) ) ) ), N3, 1 );

% fft, loop over each tensor component
C2 = zeros( N3, ndim2 );
for temp = 1 : ndim2
    C2( :, temp ) = fftfem( B2( :, temp ) );
end

% multiply Ghat4
D2 = ddot42_2( Ghat4, C2, N3 );
% D2 = ddot42_2( Ghat4, C2, N3 );
% D22 = zeros(N3,9);
% for temp = 1:N3
%     D22( temp, : ) = ddot42( Ghat4(temp,:), C2(temp,:) );
% end

% inverse fft, loop over each tensor component
G_K_dF = zeros( N3, ndim2 );
for temp = 1 : ndim2
    G_K_dF( :, temp ) = real( ifftfem( D2( :, temp ) ) );
end
% =========================================================================
% ========================== precondition =================================
% =========================================================================
% for temp = 1:ndim2
%     G_K_dF(:,temp) = G_K_dF(:,temp) ./ cg_precondition;
% end
% G_K_dF = ddot42_2(cg_precondition, G_K_dF, N3);
%
G_K_dF = ddot42_2( K_pre, G_K_dF, N3 );
G_K_dFvector = reshape( G_K_dF, ndim2 * N3, 1 );
%
% =========================================================================
% =========================================================================
% =========================================================================

end