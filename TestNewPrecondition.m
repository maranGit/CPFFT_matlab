clear
clc
%
%
%
load('Precondition_Workspace.mat');
load('K4.mat'); % load Ghat4, K4, b
G_dP = @(A2) G_convolution( A2, Ghat4, N, ndim );
G_K_dF = @(dFm) G_convolution(dFm, Ghat4, N, ndim, K4);
%
%                       original solution
%
% dFm = pcg_RM( G_K_dF, b, 1.0e-8, 1000 );
%
%               form precondition matrix (diagonal)
%
% phase = zeros(N, N, N);
% phase(3:7, 1:4, 3:7) = 1;
% phase = reshape( phase, [], 1 );
% alpha = 1.1;
% K4_bar = mean(K4);
% cg_precondition = zeros(N3,81);
% cg_precondition(1,:) = K4(1,:);
% cg_precondition(N3,:) = K4(N3,:);
% for temp = 2:(N3-1)
%     cg_precondition(temp,:) = mean(K4(temp-1:temp+1,:));
% end
% cg_precondition_inv = zeros(N3,81);
% for temp = 1:N3
%     temp2 = reshape(cg_precondition(temp,:),9,9);
%     cg_precondition_inv(temp,:) = reshape(inv(temp2),1,81);
% end
cg_precondition_inv = testa;
cg_precondition = inv(testa);
%
%                        RHS precondition
%
b_precondition = cg_precondition * b;
% b_precondition = reshape(b_precondition,343,9);
% for temp = 1:9
%     b_precondition(:,temp) = b_precondition(:,temp) ./ cg_precondition;
% end
% b_precondition = ddot42_2(cg_precondition,b_precondition, N3);
% b_precondition = reshape(b_precondition,343*9,1);
%
%                        updated solution
%
G_K_precondition = @(dFm2) G_precondition(dFm2, Ghat4, N, ndim, K4, cg_precondition, cg_precondition_inv);
dFm_precondition = pcg_RM( G_K_precondition, b_precondition, 1.0e-8, 1000 );
%
%                        recover solution
%
dFm_precondition = cg_precondition_inv * dFm_precondition;
% for temp = 1:9
%     dFm_precondition(:,temp) = dFm_precondition(:,temp) .* cg_precondition;
% end
% dFm_precondition = ddot42_2(cg_precondition_inv,dFm_precondition,N3);
%
%                            compare
%
% max(abs(dFm_precondition(:) - dFm))