function [P, A, history1] = rstgp1(Fn1, Fn, history, history1, mateprop)
% Ran Ma
% 03/19/2018
% call material subroutine
% then pull back stress and stiffness to reference configuration
%
zero = 0;
beta = mateprop(6);
PatchE = mateprop(4);
Patchv = mateprop(5);
tan_e = mateprop(7);
yld_n1 = mateprop(8);
ym_n1 = PatchE;
ym_n = PatchE;
nu_n1 = Patchv;
nu_n = Patchv;
hprime_n1 = tan_e*ym_n1/(ym_n1 - tan_e);

% for ii = 1:N3

% F to uddt
%     currF = transpose(reshape(F(ii,:),3,3));
%     oldF = transpose(reshape(F_old(ii,:),3,3));
fnh = 0.5 * (Fn1 + Fn);
dfn = Fn1 - Fn;
fnhinv = inv(fnh);
fn1inv = inv(Fn1);
detF = det(Fn1);
[rnh, ~] = poldec(fnh);
[rn1, ~] = poldec(Fn1);

ddt = dfn * fnhinv;
ddt = 0.5 * (ddt + ddt');
uddt = transpose(rnh) * ddt * rnh;
uddt_voigt = transpose(uddt([1,5,9,4,8,7]));
uddt_voigt(4:6) = uddt_voigt(4:6) * 2;

% update stress and stiffness
currhist = history(1:11);
cgn = history(12:20);

[cgn1,currhist1,rtse,yield] = ...
    mm01(ym_n1,nu_n1,beta,hprime_n1,yld_n1,cgn,uddt_voigt,currhist,ym_n,nu_n);

history1(12:20) = cgn1;
history1(1:11) = currhist1;

%% sigma to P (P = J*sigma*F^{-T})
urcs = [cgn1(1),cgn1(4),cgn1(6);cgn1(4),cgn1(2),cgn1(5);cgn1(6),cgn1(5),cgn1(3)];
sigma = rn1 * urcs * transpose(rn1);
P = sigma*transpose(fn1inv)*detF;

%% dsigma/dD to dP/dF (Eqn. 7.1.90 in Simo & Hughes)
i2v = [1,4,6; 4,2,5; 6,5,3];
i2f = [1,2,3; 4,5,6; 7,8,9];

cep = cnst1(rtse,nu_n1,ym_n1,currhist1(2),currhist1(5),beta,currhist1(1),1,1,yield);
qn1 = getrm1(rn1,2);
sigma1 = qn1*cgn1(1:6);
sigma2 = [sigma1; 0; 0; 0];
cep = ctran1(cep,qn1,sigma2,1,detF,1);

A = zeros(9,9); % dP / dF
cptau = zeros(9,9); % c_abcd + tau_ac * delta_bd

% compute c_abcd + tau_ac * delta_bd
for a=1:3
    for b=1:3
        for c=1:3
            for d=1:3
                i1 = i2f(a,b);
                i2 = i2f(c,d);
                i3 = i2v(a,b);
                i4 = i2v(c,d);
                cptau(i1,i2) = cep(i3,i4);
                if b == d
                    cptau(i1,i2) = cptau(i1,i2) + sigma(a,c)*detF;
                end
%                 if a == c
%                     cptau(i1,i2) = cptau(i1,i2) - sigma(b,d)*detF;
%                 end
            end
        end
    end
end

% compute A_aBcD = fn1inv_Bb * cptau_abcd * fn1inv_Dd
for a=1:3
    for B=1:3
        for c=1:3
            for D=1:3
                i1 = i2f(a,B);
                i2 = i2f(c,D);
                A(i1,i2) = fn1inv(B,1)*fn1inv(D,1)*cptau(i2f(a,1),i2f(c,1))...
                         + fn1inv(B,1)*fn1inv(D,2)*cptau(i2f(a,1),i2f(c,2))...
                         + fn1inv(B,1)*fn1inv(D,3)*cptau(i2f(a,1),i2f(c,3))...
                         + fn1inv(B,2)*fn1inv(D,1)*cptau(i2f(a,2),i2f(c,1))...
                         + fn1inv(B,2)*fn1inv(D,2)*cptau(i2f(a,2),i2f(c,2))...
                         + fn1inv(B,2)*fn1inv(D,3)*cptau(i2f(a,2),i2f(c,3))...
                         + fn1inv(B,3)*fn1inv(D,1)*cptau(i2f(a,3),i2f(c,1))...
                         + fn1inv(B,3)*fn1inv(D,2)*cptau(i2f(a,3),i2f(c,2))...
                         + fn1inv(B,3)*fn1inv(D,3)*cptau(i2f(a,3),i2f(c,3));
            end
        end
    end
end

% end % loop over grid point

end