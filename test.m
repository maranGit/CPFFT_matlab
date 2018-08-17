clear
clc

DbarF33=rand(3,3);
% DbarF33 = [0,0.5,0; 0.5,0,0; 0,0,0];
% DbarF33 = [-0.00506384699487628,0.100508962005208,0;0.100508962005208,-0.00506384699487628,0;0,0,0];
% DbarF33=[0,1,0;0,0,0;0,0,0];
[Pbar33, C_homo] = FFT_finite_3d(DbarF33);

C_homo_fd = zeros(9,9);
index = [1,1,1,2,2,2,3,3,3; 1,2,3,1,2,3,1,2,3];
tol = 1.0e-8;
for ii = 1:9
    DbarF33_tmp = DbarF33;
    i1 = index(1,ii);
    i2 = index(2,ii);
    DbarF33_tmp(i1,i2) = DbarF33_tmp(i1,i2) + tol;
    Pbar33_tmp = FFT_finite_3d(DbarF33_tmp);
    P_diff = ( Pbar33_tmp - Pbar33 ) / tol;
    C_homo_fd(:,ii) = P_diff';
end