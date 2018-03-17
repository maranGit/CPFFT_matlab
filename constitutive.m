
% ----------------------------- CONSTITUTIVE ------------------------------
function [P_all, K4_all, history1] = constitutive(F_all, F_old, history, history1, mateprop_all)
N3 = size( F_all, 1 );
%
% phase indicator: cubical inclusion of volume fraction (9^3) / (31^3)
% phase = zeros(N, N, N);
% phase(23:31, 1:9, 23:31) = 1; % hard phase
% phase = reshape( phase, [], 1 );
%
% frequently used identity tensor and operator
% trans2 = @(A2) A2( :, [1, 4, 7, 2, 5, 8, 3, 6, 9] ); % ij -> ji
% dot22 = @(A2, B2) reshape( reshape(B2, 3, 3) * reshape(A2, 3, 3), 1, 9 ); % ij * jk = ik
% dot24 = @(A2, B4) transpose( reshape( reshape( transpose( B4 ), 27, 3 ) * reshape( A2, 3, 3 ), 9, 9 ) ); % ij * jkmn = ikmn
% dot42 = @(A4, B2) transpose( reshape( reshape(B2, 3, 3) * reshape( transpose( A4 ), 3, 27 ), 9, 9 ) ); % ijkl * lm = ijkm
% ddot44 = @(A4, B4) A4 * B4( [ 1, 4, 7, 2, 5, 8, 3, 6, 9 ], : ); % ijkl * lkmn = ijmn
% I2 = [1, 0, 0, 0, 1, 0, 0, 0, 1]; % delta_ij
% I4 = zeros( 9, 9 ); % delta_il * delta_jk
% I4( [1, 13, 25, 29, 41, 53, 57, 69, 81] ) = 1;
% I4rt = eye( 9 ); % delta_ik * delta_jl = I4lt
% I4s = 0.5 * ( I4 + I4rt ); % ( I4 + I4rt ) / 2
% II = zeros(9, 9); % delta_ij * delta_kl
% II( [1, 5, 9, 37, 41, 45, 73, 77, 81] ) = 1;
%
% loop over integration point to update stress and consistent stiffness
P_all = zeros(N3, 9); % 1st P-K stress of each integration point
K4_all = zeros(N3, 81); % consistent stiffness of each integration point
for temp = 1 : N3 % loop over all integration point
%     if phase( temp )
%         K = 8.33;
%         mu = 3.86;
%     else
%         K = 0.833;
%         mu = 0.386;
%     end
%     F = F_all( temp, : ); % local F [ 11, 12, 13, 21, 22, 23, 31, 32, 33 ]
%     E = 0.5 * ( dot22( trans2(F), F ) - I2 );
%     C4 = K * II + 2 * mu * ( I4s - II / 3 );
%     S = transpose( C4 * transpose( E ) );
%     P = dot22( F, S );
%     K4 = dot24( S, I4 ) + ddot44( ddot44( I4rt, dot42( dot24( F, C4 ), trans2(F) ) ), I4rt );

    mateprop = mateprop_all(:,temp);
    currhist1 = history1(:,temp);
    currhist = history(:,temp);
    Fn1 = transpose(reshape(F_all(temp,:),3,3));
    Fn = transpose(reshape(F_old(temp,:),3,3));
    
    stiff1 = 0;1;
if stiff1
    % analytical tangent stiffness
    [P, A, currhist1] = fftUpdate(Fn1, Fn, currhist, currhist1, mateprop);
    history1(:,temp) = currhist1;
    
else
    
    % finite difference stiffness matrix
    A = zeros(9,9);
    tol = 1.0e-8;
    temp1 = 1;
    [P, ~, historyBackup] = fftUpdate(Fn1, Fn, currhist, currhist1, mateprop);
    for a = 1:3
        for b = 1:3
            Fn1_inc = Fn1;
            Fn1_inc(a,b) = Fn1(a,b) + tol;
            [P1, ~, ~] = fftUpdate(Fn1_inc, Fn, currhist, currhist1, mateprop);
            A(:,temp1) = reshape(transpose(P1-P)/tol,9,1);
            temp1 = temp1 + 1;
        end
    end
    history1(:,temp) = historyBackup;
end
    
    K4 = reshape(transpose(A),1,81);
    
    P_all( temp, : ) = reshape(transpose(P),1,9);
    K4_all( temp, : ) = K4;
end
end