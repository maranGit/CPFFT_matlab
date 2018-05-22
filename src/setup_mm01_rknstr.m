function [mateprop,history,history1] = setup_mm01_rknstr(N, matList, matprp)
% Ran Ma
% 03/19/2018
% initialize material parameters for mm01

N3 = N*N*N;

zero = 0;
mateprop = ones(8,N3);

for ii = 1:N3
    currmat = matList(ii);
    ym_n1 = matprp(1,currmat);
    tan_e = matprp(4,currmat);
    sigma_o = matprp(5,currmat);
    hprime = tan_e.*ym_n1./(ym_n1 - tan_e);
    root3 = sqrt(3);
    kn = sigma_o / root3;
    
    mateprop(4,ii) = matprp(1,currmat);
    mateprop(5,ii) = matprp(2,currmat);
    mateprop(6,ii) = matprp(3,currmat);
    mateprop(7,ii) = matprp(4,currmat);
    mateprop(8,ii) = matprp(5,currmat);

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