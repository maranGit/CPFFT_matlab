function Ghat4 = formG(N, ndim)

ndim2 = ndim ^ 2;
ndim3 = ndim ^ 3;
ndim4 = ndim ^ 4;
N3 = N ^ 3;

% identity tensor
delta = eye( ndim );

% projection operator
q = zeros( N3, ndim ); % frequency vector [ N1, N2, N3 ]
for ii = 0 : ( ndim - 1 )
    q( :, ii + 1 ) = repmat( reshape( repmat( -( N-1 ) / 2 : ( N-1 ) / 2, ...
    N^( ndim - ii - 1 ), 1 ), [], 1 ), N^ii, 1 );
end
q_dot_q = sum( q .* q, 2 );
q_dot_q( abs( q_dot_q ) < 1e-3 ) = Inf; % zero frequency -> mean
Ghat4 = zeros( N3, ndim4 );
for temp = 0 : ( ndim4 - 1 ) % loop over all integration point
    ii = mod( fix( temp / ndim3 ), ndim ) + 1; % i of Ghat4_ijlm
    jj = mod( fix( temp / ndim2 ), ndim ) + 1; % j of Ghat4_ijlm
    ll = mod( fix( temp / ndim ), ndim ) + 1; % l of Ghat4_ijlm
    mm = mod( temp, 3 ) + 1; % m of Ghat4_ijlm
    Ghat4( :, temp + 1 ) = delta(ii, ll) * q( :, jj ) .* q( :, mm ) ./ q_dot_q;
end

end