% multi elements test
ndm = 2;
nen = 4;
numnp = 16;
numel = 9;
ndf = 2;
Coordinates = [0, 0; 1, 0; 2, 0; 3, 0;
    0, 1; 1, 1; 2, 1; 3, 1;
    0, 2; 1, 2; 2, 2; 3, 2;
    0, 3; 1, 3; 2, 3; 3, 3];
NodesOnElement = [0 1 5 4;
              1 2 6 5;
              2 3 7 6;
              4 5 9 8;
              5 6 10 9;
              6 7 11 10;
              8 9 13 12;
              9 10 14 13;
              10 11 15 14];
NodesOnElement = NodesOnElement + 1;
RegionOnElement = ones(9, 1);
% eType = 'ele3';
FaceBC = 0;
NodeBC = [1, 1, 0;
       1, 2, 0;
       5, 1, 0;
       9, 1, 0;
       13, 1, 0;
       2, 2, 0;
       3, 2, 0;
       4, 2, 0;
       4, 1, 0.6;
       8, 1, 0.6;
       12, 1, 0.6;
       16, 1, 0.6];
numBC = 12;
NBC = [13, 2, 0;
       14, 2, 0;
       15, 2, 0;
       16, 2, 0];
nummat = 1;
mults = 0.02:0.02:1;
% mat = {'VonMises_finite',[12000;0.3;100;1000;3]};

nonlin = 1;
ProbType = [numnp, numel, nummat, ndm, ndf, nen];
AlgoType = [0; 1; 0];
MatTypeTable = [1; 46; 1];
OptFlag = [0 1 1 0 0 1 1 0]';
MateT = [1, 1, 1, 12000, 0.3, 0.5, 1000, 100];
%
stepmax = 50;
datastep = 50; % size of output arrays, can be larger than stepmax
 
SHList = 1:datastep; % Ask for stress at all time steps
SEHList = 1:datastep; % Ask for stress at all time steps
 
%%%%% FIX this strain rate and tstep
% strainstep = diplac*s_del_a; % strain applied per step
% rate = 1.2e-7; % desired strain rate
tstep = 1; % makes the rate be 'rate' s-1
 
itermax = 8;
Residratio = 10^-7;
reststep = 200; % dump data at steps equalting multiples of this #
LieHist = 0;
ReMeth = 0;
 
restartmat = 1;
printRnorm = 1;
hrstore = 1;
 
numsubcycles = 0;