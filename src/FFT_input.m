

nstep = 10;
straininc = 0.001;1;
N = 7; % number of voxels (assumed equal for all directions)
matprp = zeros(300,500);
phase = zeros(N, N, N);
temp = round(N/2);
phase(N-temp:N, 1:temp, N-temp:N) = 1; % hard phase
matList = reshape( phase, [], 1 ) + 1;
assigned = zeros(500,1);
assigned(1:2) = 1;
matprp(1:5,1) = [12000;0.3;0.5;1000;100];
matprp(1:5,2) = [24000;0.3;0.5;1000;200];