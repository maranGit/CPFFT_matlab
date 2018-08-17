function [K_pre, K_pre_inv] = formPrecondition(K, isw, alpha)

lastwarn('');
N = size(K,1);
M = size(K,2);
if M ~= 81
    error('>>> Error: K should have 81 components');
end

switch isw
    case 1
        K_pre = Kavm(K, alpha);
    case 2
        K_pre = Kavg(K, alpha);
    case 3
        K_pre = Kavh(K, alpha);
    otherwise
        error('>>> Error: Unknow option');
end

% inverse of precondition matrix
K_pre_inv = zeros(N,81);
for temp = 1:N
    K_temp = reshape(K_pre(temp,:),9,9);
    K_temp_inv = pinv(K_temp);

    % check if K_temp_inv is singular
    [~, msgid] = lastwarn;
    if strcmp(msgid,'MATLAB:singularMatrix') || strcmp(msgid,'MATLAB:nearlySingularMatrix')
        error('>>> Error: precondition matrix is singular');
    end

    K_pre_inv(temp,:) = reshape(K_temp_inv,1,81);
end
end

%% algorithmic average
function Kava = Kavm(K, alpha)
N = size(K,1);
K_bar = repmat(mean(K),N,1);
Kava = K_bar + (K - K_bar) * alpha;
end

%% geometric average
function Kavg2 = Kavg(K, ~)
N = size(K,1);
K_bar = repmat(mean(K),N,1);
Kboth = [reshape(K, 81*N, 1), reshape(K_bar, 81*N, 1)];
Kavg = reshape(geomean(Kboth,2), N, 81);
Kavg2 = Kavg + repmat(mean(K-Kavg), N, 1);
end

%% harmonic average
function Kavh2 = Kavh(K, ~)
N = size(K,1);
K_bar = repmat(mean(K),N,1);
Kboth = [reshape(K, 81*N, 1), reshape(K_bar, 81*N, 1)];
Kavh = reshape(harmmean(Kboth,2), N, 81);
Kavh2 = Kavh + repmat(mean(K-Kavh), N, 1);
end