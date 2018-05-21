% Ran Ma
% 12/19/2017
% 
% De Geus, T. W. J., et al. "Finite strain FFT-based non-linear solvers 
% made simple." Computer Methods in Applied Mechanics and Engineering 
% 318 (2017): 412-430.
% 
% column arrangement
% 2nd order tensor: 11, 12, 13, 21, 22, 23, 31, 32, 33
% 4th order tensor: 1111, 1112, 1113, ..., 1133, 1211, ..., 3332, 3333
%
% row arrangement ( N1, N2, N3 )
% (111), (112), ..., (11N), (121), (122), ..., (NN1), (NN2), ..., (NNN)
% 
% update 12/21/2017: verified against hyperelasticity python code
%                    both stress and residuals are exactly the same
%
clear
clc

FFT_input

FFT_init

FFT_nr3

% FFT_post