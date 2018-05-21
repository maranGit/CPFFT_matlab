function [mateprop,history,history1] = fftInit(N, matList, matprp)
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
mateprop(4,temp2) = 24000.0;
mateprop(5,temp2) = 0.3;
mateprop(6,temp2) = 0.5;
mateprop(7,temp2) = 1000.0;
mateprop(8,temp2) = 200.0;

zero = 0;
for ii = 1:N3
    currmat = matList(ii);
    ym_n1 = matprp(1,currmat);
    tan_e = matprp(4,currmat);
    sigma_o = matprp(5,currmat);
    hprime = tan_e.*ym_n1./(ym_n1 - tan_e);
    root3 = sqrt(3);
    kn = sigma_o / root3;

    % initialize history array
    history = zeros(20,N3);
    history(1,ii) = zero;
    history(2,ii) = kn;
    history(3,ii) = zero;
    history(4,ii) = 3.0;
    history(5,ii) = hprime;
    history(6:20,ii) = zero;
    
end
history1 = history;

end