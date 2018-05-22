function [P_all, K4_all, history1] = drive_sig_up(F, F_old, history, history1, mateprop)
% Ran Ma
% 03/19/2018
% take P and A from lower level function
% then 'assemble' P and A to P_all and K4_all
%
N3 = size( F, 1 );
%
% loop over integration point to update stress and consistent stiffness
P_all = zeros(N3, 9); % 1st P-K stress of each integration point
K4_all = zeros(N3, 81); % consistent stiffness of each integration point

for temp = 1 : N3 % loop over all integration point

    mateprop_local = mateprop(:,temp);
    currhist1 = history1(:,temp);
    currhist = history(:,temp);
    Fn1_local = transpose(reshape(F(temp,:),3,3));
    Fn_local = transpose(reshape(F_old(temp,:),3,3));
    
    stiff_flag = 2;0;1;
    
    switch stiff_flag
        case 0
            % analytical tangent stiffness
            [P_local, A_local, currhist1] = rstgp1(Fn1_local, Fn_local, currhist, currhist1, mateprop_local);
            history1(:,temp) = currhist1;
        case 1
            % tangent stiffness from finite difference method
            A_local = zeros(9,9);
            tol = 1.0e-8;
            temp1 = 1;
            [P_local, ~, historyBackup] = rstgp1(Fn1_local, Fn_local, currhist, currhist1, mateprop_local);
            for a = 1:3
                for b = 1:3
                    Fn1_inc = Fn1_local;
                    Fn1_inc(a,b) = Fn1_local(a,b) + tol;
                    [P1, ~, ~] = rstgp1(Fn1_inc, Fn_local, currhist, currhist1, mateprop_local);
                    A_local(:,temp1) = reshape(transpose(P1-P_local)/tol,9,1);
                    temp1 = temp1 + 1;
                end
            end
            history1(:,temp) = historyBackup;
            1;
        case 2
            % verify analytical tangent stiffness against finite difference
            % (1) finite difference
            A_local = zeros(9,9);
            tol = 1.0e-8;
            temp1 = 1;
            [P_local, ~, historyBackup] = rstgp1(Fn1_local, Fn_local, currhist, currhist1, mateprop_local);
            for a = 1:3
                for b = 1:3
                    Fn1_inc = Fn1_local;
                    Fn1_inc(a,b) = Fn1_local(a,b) + tol;
                    [P1, ~, ~] = rstgp1(Fn1_inc, Fn_local, currhist, currhist1, mateprop_local);
                    A_local(:,temp1) = reshape(transpose(P1-P_local)/tol,9,1);
                    temp1 = temp1 + 1;
                end
            end
            
            % (2) analytical
            [~, A2, ~] = rstgp1(Fn1_local, Fn_local, currhist, currhist1, mateprop_local);
            
            % (3) compare difference
%             A_diff = ( abs(A - A2) ./ ( abs(A) + abs(A2) ) );
            A_diff = sumsqr(A_local - A2) / sumsqr(A_local + A2);
%             max_A_diff = max(A_diff(:),[],'omitnan');
            if A_diff > 0.01
                error('>>> Error: analytical stiffness matrix is not accurate enough');
            end
            
            % update history
            history1(:,temp) = historyBackup;
            
        otherwise
            error('>>> Error: Unrecogonized option for stiff_flag');
    end
    
    K4 = reshape(transpose(A_local),1,81);
    
    P_all( temp, : ) = reshape(transpose(P_local),1,9);
    K4_all( temp, : ) = K4;
end
end