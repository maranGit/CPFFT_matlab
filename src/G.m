
% ------------------------- PROJECTION OPERATION --------------------------
function G_K_dFvector = G(A2vector, Ghat4, N, ndim, K4)
%
% problem size
ndim2 = ndim ^ 2;
N3 = N ^ 3;

% 
A2 = reshape( A2vector, N3, ndim2 );

% tensor operations/products
trans2 = @(A2) A2( :, [1, 4, 7, 2, 5, 8, 3, 6, 9] ); % ij -> ji
ddot42 = @(A4, B2) transpose( reshape( sum( reshape( transpose( A4 .* ...
     repmat( trans2(B2), 1, ndim2 ) ), ndim2, ndim2 * N3 ) ), ndim2, N3) );
% C_ij = A4_ijkl * B2_lk

if nargin == 4 % A2 is 1st P-K stress, don't multiply K4
    B2 = A2;
elseif nargin == 5 % A2 is strain, multiply K4
    B2 = trans2( ddot42( K4, trans2( A2 ) ) );
end

% (inverse) Fourier transform (for each tensor component in each direction)
fftfem = @(x) reshape( fftshift( fftn( ifftshift( reshape( x, N, N, N ) ) ) ), N3, 1 );
ifftfem = @(x) reshape( fftshift( ifftn( ifftshift( reshape( x, N, N, N ) ) ) ), N3, 1 );

% fft, loop over each tensor component
C2 = zeros( N3, ndim2 );
for temp = 1 : ndim2
    C2( :, temp ) = fftfem( B2( :, temp ) );
end

% multiply Ghat4
D2 = ddot42( Ghat4, C2 );

% inverse fft, loop over each tensor component
G_K_dF = zeros( N3, ndim2 );
for temp = 1 : ndim2
    G_K_dF( :, temp ) = real( ifftfem( D2( :, temp ) ) );
end

% reshape to vector form ( ndim2 x N3, 1)
G_K_dFvector = reshape( G_K_dF, ndim2 * N3, 1 );
end