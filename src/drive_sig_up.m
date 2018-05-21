function [P_all, K4_all, history1] = drive_sig_up(F_all, F_old, history, history1, mateprop_all)
% Ran Ma
% 03/19/2018
% take P and A from lower level function
% then 'assemble' P and A to P_all and K4_all
%
N3 = size( F_all, 1 );
%
% loop over integration point to update stress and consistent stiffness
P_all = zeros(N3, 9); % 1st P-K stress of each integration point
K4_all = zeros(N3, 81); % consistent stiffness of each integration point

for temp = 1 : N3 % loop over all integration point

    mateprop = mateprop_all(:,temp);
    currhist1 = history1(:,temp);
    currhist = history(:,temp);
    Fn1 = transpose(reshape(F_all(temp,:),3,3));
    Fn = transpose(reshape(F_old(temp,:),3,3));
    
    stiff_flag = 2;0;1;
    
    switch stiff_flag
        case 0
            % analytical tangent stiffness
            [P, A, currhist1] = rstgp1(Fn1, Fn, currhist, currhist1, mateprop);
            history1(:,temp) = currhist1;
        case 1
            % tangent stiffness from finite difference method
            A = zeros(9,9);
            tol = 1.0e-8;
            temp1 = 1;
            [P, ~, historyBackup] = rstgp1(Fn1, Fn, currhist, currhist1, mateprop);
            for a = 1:3
                for b = 1:3
                    Fn1_inc = Fn1;
                    Fn1_inc(a,b) = Fn1(a,b) + tol;
                    [P1, ~, ~] = rstgp1(Fn1_inc, Fn, currhist, currhist1, mateprop);
                    A(:,temp1) = reshape(transpose(P1-P)/tol,9,1);
                    temp1 = temp1 + 1;
                end
            end
            history1(:,temp) = historyBackup;
            1;
        case 2
            % verify analytical tangent stiffness against finite difference
            % (1) finite difference
            A = zeros(9,9);
            tol = 1.0e-8;
            temp1 = 1;
            [P, ~, historyBackup] = rstgp1(Fn1, Fn, currhist, currhist1, mateprop);
            for a = 1:3
                for b = 1:3
                    Fn1_inc = Fn1;
                    Fn1_inc(a,b) = Fn1(a,b) + tol;
                    [P1, ~, ~] = rstgp1(Fn1_inc, Fn, currhist, currhist1, mateprop);
                    A(:,temp1) = reshape(transpose(P1-P)/tol,9,1);
                    temp1 = temp1 + 1;
                end
            end
            
            % (2) analytical
            [~, A2, ~] = rstgp1(Fn1, Fn, currhist, currhist1, mateprop);
            
            % (3) compare difference
%             A_diff = ( abs(A - A2) ./ ( abs(A) + abs(A2) ) );
            A_diff = sumsqr(A - A2) / sumsqr(A + A2);
%             max_A_diff = max(A_diff(:),[],'omitnan');
            if A_diff > 0.01
                error('>>> Error: analytical stiffness matrix is not accurate enough');
            end
            
            % update history
            history1(:,temp) = historyBackup;
            
        otherwise
            error('>>> Error: Unrecogonized option for stiff_flag');
    end
    
    K4 = reshape(transpose(A),1,81);
    
    P_all( temp, : ) = reshape(transpose(P),1,9);
    K4_all( temp, : ) = K4;
end
end