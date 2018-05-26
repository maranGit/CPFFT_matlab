%
%      last modified:  2/17/2017 rhd
%
%            Ran: number of gauss point in each fft 'element'
fftngp=1;

mxnod=3000000;
mxel=3000000;
mxndel=20;
max_threads=128;
mxlbel=15;
mxndof=3;
mxgp=14;
nstr=6;
mxstep=100000;
mxlc=5000;
mxmat=500;
mxlsz=mxnod/3;
mxelmp=30;
mxstmp=30;
mxstc=10;
mxoupr=30;
ntrc=10;
mxconn=100;
mxvl=128;
mxnmgp=10;
mxblsz=128;
mxnmbl=20000;
nparam=3;
mxelpr=42;
vclim=28;
mxmtpr=300;
mxsepr=1024;
two16=65536;
ndim=3;
max_tables=20;
mxtim=50;
nstrs=9;
max_procs=2048;
mxndpr=8;
max_crystals=1000;
max_slip_sys=48;
max_uhard=48;
max_surfaces=100;
max_tied_sets=100;
max_packet_types=200;
max_user_lists=100;
mxcvtests=10; nonlocal_shared_state_size=12;
max_interface_props=80;

mxdof=mxnod*mxndof;
trsz=mxnod*mxndof;
mxncor=3*mxnod;
mxecor=3*mxndel;
mxedof=mxndof*mxndel;
mxtgp=mxel*mxgp;
mxndldcm=mxndof+1;
mxtnsz=ndim*ndim;mxndlm=(mxnod/31)+1;
mxutsz=(mxedof*mxedof+mxedof)/2;
mxnusz=(mxndel*mxndel+mxndel)/2;
max_mpc=mxnod/10;

mxoupt=mxndel;