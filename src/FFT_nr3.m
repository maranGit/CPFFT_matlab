
% function for the projection 'G', and the proeuct 'G : K^LT : (delta F)^T'
G_dP = @(A2) G( A2, Ghat4, N, ndim ); % A2 can be either vector or matrix

% --------------------------- NEWTON ITERATIONS ---------------------------
% initialize deformation gradient, and stress/stiffness
F = zeros( N3, 3 * ndim );
F( :, [1, 5, 9] ) = 1.0; % 11, 12, 13, 21, 22, 23, 31, 32, 33
F_old = F;

% initial internal variables and stiffness matrix
[mateprop,history,history1] = fftInit(N, matList, matprp);
[~, K4, ~] = drive_sig_up(F, F_old, history, history1, mateprop);
[mateprop,history,history1] = fftInit(N, matList, matprp);

% set macroscopic loading
barF = zeros( N3, 3 * ndim );
barF( :, [1, 5, 9] ) = 1.0; % 11, 12, 13, 21, 22, 23, 31, 32, 33
barF_old = barF;
stressStrain = zeros(nstep,2);

for step = 1:nstep
    
    fprintf('>>> Now starting step %i \n', step)
    iiter = 0;
    
    % set macroscopic deformation gradient (pure-shear)
    barF = zeros( N3, 3 * ndim );
    barF( :, [1, 5, 9] ) = 1.0; % 11, 12, 13, 21, 22, 23, 31, 32, 33
    barF(:, 1) = 1 + step * straininc;
    barF(:, 5) = 1.0 / (1.0 + step * straininc);
    
    % store normalization
    Fn = sqrt( sumsqr( F ) );
    
    % first iteration residual; distribute "barF" over grid using "K4"
    G_K_dF = @(dFm) G(dFm, Ghat4, N, ndim, K4);
    b = - G_K_dF( barF - barF_old );
    F = F + barF - barF_old;
    
    % iterate as long as the iterative update does not vanish
    while 1
        fprintf('>>> Now starting iteration %i \n', iiter)
        dFm = pcg( G_K_dF, b, 1.0e-10, N3 ); % solve linear system using CG
        F = F + reshape( dFm, N3, ndim2 ); % update DOFs
        [P, K4, history1] = drive_sig_up(F, F_old, history, history1, mateprop); % new residual stress and tangent
        b = - G_dP( P ); % convert residual stress to residual
        G_K_dF = @(dFm) G(dFm, Ghat4, N, ndim, K4); % update K4 in anonymous function
        res = sqrt( transpose(dFm) * dFm );
        
        fprintf('relative residual is %6.3e \n', res / Fn) % print residual
        if res/Fn < 1.0e-5 && iiter > 0
            break;
        end
        if iiter == 10
            error('>>>Error: Newton loop does not converge.');
        end
        iiter = iiter + 1;
    end
    
    stressStrain(step,:) = [F(1,1),P(1,1)];
    
%     fprintf('hardening variable is %7.3e \n',history1(1))
    
    if iiter == 10
        fprintf('>>> Error: do not converge')
        break;
    end
    
    % update current state
    F_old = F;
    barF_old = barF;
    history = history1;
    
end