
% --------------------------------- GRID ----------------------------------
ndim = 3; % number of dimensions, should be fixed as 3 for now
ndim2 = ndim * ndim;
ndim3 = ndim2 * ndim;
ndim4 = ndim3 * ndim;
N3 = N ^ 3;

Ghat4 = formG(ndim,N,N3);

[nelblk, elblks] = inelbk( matList, N3 );

urcs_blk_list = zeros(nelblk,1);
eps_blk_list = zeros(nelblk,1);
rot_blk_list = zeros(nelblk,1);
history_blk_list = zeros(nelblk,1);
cep_blk_list = zeros(nelblk,1);

urcs_n_blocks(1:nelblk) = blocks_ptr_type();
urcs_n1_blocks(1:nelblk) = blocks_ptr_type();
eps_n_blocks(1:nelblk) = blocks_ptr_type();
eps_n1_blocks(1:nelblk) = blocks_ptr_type();
rot_n_blocks(1:nelblk) = blocks_ptr_type();
rot_n1_blocks(1:nelblk) = blocks_ptr_type();
rot_n_blocks(1:nelblk) = blocks_ptr_type();
rot_n1_blocks(1:nelblk) = blocks_ptr_type();
history_blocks(1:nelblk) = blocks_ptr_type();
history_n1_bocks(1:nelblk) = blocks_ptr_type();
cep_blocks(1:nelblk) = blocks_ptr_type();