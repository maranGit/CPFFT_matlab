function Pbar33 = FFT_finite_3d(DbarF33, isprecondition)
% Ran Ma
% 12/19/2017
% 
% De Geus, T. W. J., et al. "Finite strain FFT-based non-linear solvers 
% made simple." Computer Methods in Applied Mechanics and Engineering 
% 318 (2017): 412-430.
% 
% column arrangement
% 2nd order tensor: 11, 12, 13, 21, 22, 23, 31, 32, 33
% 4th order tensor: 1111, 1112, 1113, ..., 1133, 1211, ..., 3332, 3333
%
% row arrangement ( N1, N2, N3 )
% (111), (112), ..., (11N), (121), (122), ..., (NN1), (NN2), ..., (NNN)
% 
% update 12/21/2017: verified against hyperelasticity python code
%                    both stress and residuals are exactly the same
%
% clear
% clc
% --------------------------------- GRID ----------------------------------
ndim = 3; % number of dimensions, should be fixed as 3 for now
ndim2 = ndim ^ 2;
N = 7; % number of voxels (assumed equal for all directions)
N3 = N ^ 3;

% -------------------------- PRECONDITION OPTION --------------------------

poption = 1;
alpha = 0.9999;

% -------------------- PROJECTION, TENSORS, OPERATIONS --------------------
% tensor operations/products
% trans2 = @(A2) A2( :, [1, 4, 7, 2, 5, 8, 3, 6, 9] ); % transfer operation of 2nd order tensor

Ghat4 = formG(N, ndim);
% Ghat4 = formG(N, ndim);

% function for the projection 'G', and the product 'G : K^LT : (delta F)^T'
% G_dP = @(A2) G_convolution( A2, Ghat4, N, ndim ); % A2 can be either vector or matrix

% --------------------------- NEWTON ITERATIONS ---------------------------
% initialize deformation gradient, and stress/stiffness
F = zeros( N3, 3 * ndim );
F( :, [1, 5, 9] ) = 1.0; % 11, 12, 13, 21, 22, 23, 31, 32, 33
[~, K4] = constitutive(F);

nsteps = 1;
for step = 1:nsteps
% set macroscopic loading
% DbarF = zeros( N3, 3 * ndim );
% DbarF(:, 2) = 1.0/nsteps;
DbarF = repmat(reshape(DbarF33'/nsteps,1,9),N3,1);

% initial residual: distribute "barF" over grid using K4
% G_K_dF = @(dFm) G_convolution(dFm, Ghat4, N, ndim, K4);
%=========================================================================
%                            precondition
%=========================================================================
if isprecondition == 0
    K_pre = repmat(reshape(eye(9),1,81),N3,1); K_pre_inv = K_pre; % no precondition
elseif isprecondition == 1
    [K_pre,K_pre_inv] = formPrecondition(K4, poption, alpha); % precondition
elseif isprecondition == 3
    K_pre = repmat(reshape(eye(9),1,81),N3,1); K_pre_inv = K_pre; % no precondition
    [K4, ~] = formPrecondition(K4, poption, alpha); % modified Newton iteration
end
G_K_dF = @(dFm) G_precondition(dFm, Ghat4, N, ndim, K4, K_pre, K_pre_inv);
b = -G_convolution(DbarF, Ghat4, N, ndim, K4);
b = ddot42_2( K_pre, b, N3 );
%=========================================================================
%=========================================================================
%=========================================================================
F = F + DbarF;
Fn = sqrt( sumsqr( F ) );
res = Fn;
iiter = 0;

% iterate as long as the iterative update does not vanish
while ( res / Fn > 1.0e-5 )
    dFm = pcg_RM( G_K_dF, b, 1.0e-8, 1000, K_pre_inv ); % solve linear system using CG
    dFm = ddot42_2( K_pre_inv, dFm, N3 );
    F = F + reshape( dFm, N3, ndim2 ); % update DOFs
    [P, K4] = constitutive(F); % new residual stress and tangent
% =========================================================================
%                            precondition
% =========================================================================
    if isprecondition == 0
        K_pre = repmat(reshape(eye(9),1,81),N3,1); K_pre_inv = K_pre; % no precondition
    elseif isprecondition == 1
        [K_pre,K_pre_inv] = formPrecondition(K4, poption, alpha); % precondition
    elseif isprecondition == 3
        K_pre = repmat(reshape(eye(9),1,81),N3,1); K_pre_inv = K_pre; % no precondition
        [K4, ~] = formPrecondition(K4, poption, alpha); % modified Newton iteration
    end
    b = -G_convolution(P, Ghat4, N, ndim);
    b = ddot42_2( K_pre, b, N3 );
    G_K_dF = @(dFm) G_precondition(dFm, Ghat4, N, ndim, K4, K_pre, K_pre_inv);
% =========================================================================
% =========================================================================
% =========================================================================
    res = sqrt( transpose(dFm) * dFm );
    sprintf('relative residual is %10.9e', res / Fn) % print residual
    iiter = iiter + 1;
end % N-R loop
end % end all steps

% get homogenized first PK stress
Pbar33 = mean(P);

% get homogenized tangent stiffness
% C_homo = tangent_homo(Ghat4, K4, N);
% C_homo = 0;

end
%% ----------------PROBLEM DEFINITION / CONSTITUTIVE MODEL-----------------
function [P_all, K4_all] = constitutive(F_all) % verified against python
N3 = size( F_all, 1 );
N = nthroot( N3, 3 );
%
% phase indicator: cubical inclusion of volume fraction (9^3) / (31^3)
phase = zeros(N, N, N);
% phase(16:31, 1:16, 16:31) = 1; % hard phase
phase(3:7, 1:4, 3:7) = 1;
phase = reshape( phase, [], 1 );
%
% frequently used identity tensor and operator
trans2 = @(A2) A2( :, [1, 4, 7, 2, 5, 8, 3, 6, 9] ); % ij -> ji
dot22 = @(A2, B2) reshape( reshape(B2, 3, 3) * reshape(A2, 3, 3), 1, 9 ); % ij * jk = ik
dot24 = @(A2, B4) transpose( reshape( reshape( transpose( B4 ), 27, 3 ) * reshape( A2, 3, 3 ), 9, 9 ) ); % ij * jkmn = ikmn
dot42 = @(A4, B2) transpose( reshape( reshape(B2, 3, 3) * reshape( transpose( A4 ), 3, 27 ), 9, 9 ) ); % ijkl * lm = ijkm
ddot44 = @(A4, B4) A4 * B4( [ 1, 4, 7, 2, 5, 8, 3, 6, 9 ], : ); % ijkl * lkmn = ijmn
I2 = [1, 0, 0, 0, 1, 0, 0, 0, 1]; % delta_ij
I4 = zeros( 9, 9 ); % delta_il * delta_jk
I4( [1, 13, 25, 29, 41, 53, 57, 69, 81] ) = 1;
I4rt = eye( 9 ); % delta_ik * delta_jl = I4lt
I4s = 0.5 * ( I4 + I4rt ); % ( I4 + I4rt ) / 2
II = zeros(9, 9); % delta_ij * delta_kl
II( [1, 5, 9, 37, 41, 45, 73, 77, 81] ) = 1;
%
% loop over integration point to update stress and consistent stiffness
P_all = zeros(N3, 9); % 1st P-K stress of each integration point
K4_all = zeros(N3, 81); % consistent stiffness of each integration point
for temp = 1 : N3 % loop over all integration point
    if phase( temp )
        K = 0.833e-3;
        mu = 0.833e-3;
    else
        K = 0.833;
        mu = 0.386;
    end
    F = F_all( temp, : ); % local F [ 11, 12, 13, 21, 22, 23, 31, 32, 33 ]
    E = 0.5 * ( dot22( trans2(F), F ) - I2 );
    C4 = K * II + 2 * mu * ( I4s - II / 3 );
    S = transpose( C4 * transpose( E ) );
    P = dot22( F, S );
    K4 = dot24( S, I4 ) + ddot44( ddot44( I4rt, dot42( dot24( F, C4 ), trans2(F) ) ), I4rt );
    K4 = K4([1,4,7,2,5,8,3,6,9],:);
    P_all( temp, : ) = P;
    K4_all( temp, : ) = reshape( transpose(K4), 1, 81);
end
end

