function Ghat4 = formG2(N, ndim)

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

% q = q / N3;
q = q * pi / N;

k = zeros(N3,3);
temp = 0.25 * 1i * (1+exp(1i*q(:,1))) .* (1+exp(1i*q(:,2))) .* (1+exp(1i*q(:,3)));
k(:,1) = temp .* tan(q(:,1)/2);
k(:,2) = temp .* tan(q(:,2)/2);
k(:,3) = temp .* tan(q(:,3)/2);
% k = 1i*sin(q);
k_dot_k = sum(k.*conj(k),2);
k_dot_k(abs(k_dot_k)<1e-10) = inf;

Ghat4 = zeros( N3, ndim4 );
for temp = 0 : ( ndim4 - 1 ) % loop over all integration point
    ii = mod( fix( temp / ndim3 ), ndim ) + 1; % i of Ghat4_ijlm
    jj = mod( fix( temp / ndim2 ), ndim ) + 1; % j of Ghat4_ijlm
    ll = mod( fix( temp / ndim ), ndim ) + 1; % l of Ghat4_ijlm
    mm = mod( temp, 3 ) + 1; % m of Ghat4_ijlm
    Ghat4( :, temp + 1 ) = delta(ii, ll) * k( :, jj ) .* conj(k( :, mm )) ./ k_dot_k;
end

end

% 0.0266395876449033;
% -0.00204753129476961;
% 0.00447455448455661;
% 0.00847639202922296;
% 0.0507136674141935;
% 0.0405914412828459;
% 0.0721507122887331;
% 0.0352312784058994;
% -0.00757167066154966;
% 0.00650963644464527;