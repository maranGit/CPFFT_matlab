
% --------------------------------- GRID ----------------------------------
ndim = 3; % number of dimensions, should be fixed as 3 for now
ndim2 = ndim * ndim;
ndim3 = ndim2 * ndim;
ndim4 = ndim3 * ndim;
N3 = N ^ 3;

Ghat4 = formG(ndim,N,N3);