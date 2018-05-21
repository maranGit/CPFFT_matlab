clear;clc;

% problem dimension
N3 = 3;

% construct F
F_old = zeros(N3,9);
for ii = 1:N3
    F_old(:,1) = 1;
    F_old(:,2) = 0;
    F_old(:,3) = 0;
    F_old(:,4) = 0;
    F_old(:,5) = 1;
    F_old(:,6) = 0;
    F_old(:,7) = 0;
    F_old(:,8) = 0;
    F_old(:,9) = 1;
end
F = F_old;
F(1,:) = F(1,:) + rand(1,9)/10000;
F(2,:) = F(1,:);
F(3,:) = F(1,:);

% initialize material parameters and history
[mateprop,history,history1] = fftInit(N3);

% update stress and stiffness
for ii = 1:N3
    currhist1 = history1(:,ii);
    currhist = history(:,ii);
    Fn1 = transpose(reshape(F(ii,:),3,3));
    Fn = transpose(reshape(F_old(ii,:),3,3));
    [P, A, currhist1] = fftUpdate(Fn1, Fn, currhist, currhist1, mateprop);
    history1(:,ii) = currhist1;
end

%% compare with finite difference

A_fd = zeros(9,9);

tol = 1.0e-10;
temp = 1;
for a = 1:3
    for b = 1:3
        [mateprop,history,history1] = fftInit(N3);
        currhist1 = history1(:,1);
        currhist = history(:,1);
        Fn1 = transpose(reshape(F(1,:),3,3));
        Fn = transpose(reshape(F_old(1,:),3,3));
        Fn1(a,b) = Fn1(a,b) + tol;
        [P1, ~, ~] = fftUpdate(Fn1, Fn, currhist, currhist1, mateprop);
        A_fd(:,temp) = reshape(transpose(P1-P)/tol,9,1);
        temp = temp + 1;
    end
end

% A_fd = (A_fd+A_fd')/2;
A_diff = (A-A_fd)./A_fd;
max(max(A_diff))