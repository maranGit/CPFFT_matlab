function [mateprop,history,history1] = fftInit(N)
% Ran Ma
% 03/19/2018
% initialize material parameters for mm01

N3 = N*N*N;
temp = round(N/2);
phase = zeros(N, N, N);
phase(N-temp:N, 1:temp, N-temp:N) = 1; % hard phase
phase = reshape( phase, [], 1 );

% hard code material parameters
temp = [1;1;1;12000;0.3;0.5;1000;100];
mateprop = repmat(temp,1,N3);

% hard phase
temp2 = find(phase);
mateprop(8,temp2) = mateprop(8,temp2) * 2;

zero = 0;
PatchE = mateprop(4,:);
sigma_o = mateprop(8,:);
tan_e = mateprop(7,:);
ym_n1 = PatchE;
hprime = tan_e.*ym_n1./(ym_n1 - tan_e);
root3 = sqrt(3);
kn = sigma_o / root3;
Patchv = mateprop(5,:);
% nu_n1 = Patchv;
% ym_n = PatchE;
% nu_n = Patchv;

% initialize history array
history = zeros(20,N3);
history(1,:) = zero;
history(2,:) = kn;
history(3,:) = zero;
history(4,:) = 3.0;
history(5,:) = hprime;
history(6:20,:) = zero;

history1 = history;

end