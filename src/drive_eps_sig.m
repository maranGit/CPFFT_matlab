%     ****************************************************************
%     *                                                              *
%     *                  subroutine  drive_eps_sig                   *
%     *                                                              *
%     *                       written by : RM                        *
%     *                                                              *
%     *                   last modified : 4/5/2018 RM                *
%     *                                                              *
%     *      recovers all the strains, stresses                      *
%     *      for all the elements in the structure at state (n+1).   *
%     *                                                              *
%     *      Blocks are processed in parallel; using threads (OMP)   *
%     *                                                              *
%     ****************************************************************
%
function drive_eps_sig( step, iiter )
%     implicit none
%     include 'common.main'

%                         global
%     integer :: step, iiter

%                         local
%     integer :: blk, now_thread
%     logical :: debug
%     integer, external :: omp_get_thread_num
%
%     debug = .false.
%
%     call omp_set_dynamic( .false. )
%$OMP PARALLEL DO ORDERED
%$OMP&         PRIVATE( blk, now_thread )
%$OMP&         SHARED( nelblk, iiter, step )
for blk = 1:nelblk
    do_nleps_block( blk, iiter, step );
end
%$OMP END PARALLEL DO

end
%
%     ****************************************************************
%     *                                                              *
%     *                      subroutine do_nleps_block               *
%     *                                                              *
%     *                       written by : RM                        *
%     *                                                              *
%     *                   last modified : 02/26/2018 RM              *
%     *                                                              *
%     ****************************************************************
%
function do_nleps_block( blk, iter, step )
%     use fft, only: Fn1, Pn1, Fn, K4, matList, tstep, matprp
%     use elem_block_data
%     implicit none
%     include 'common.main'
%     include 'include_sig_up'
%
%             global
%     integer, intent(in) :: blk, iter, step
%
%             locals
%
%     integer :: ii, span, felem, currElem, currmat
%     integer :: info_vector(4)
%     logical :: local_debug, include_qbar
%
%     integer :: ngp, gpn
%     real(8), parameter :: zero = 0.0D0, one = 1.0D0, half = 0.5D0
%     double precision :: gp_dtemps(mxvl)
%     integer :: alloc_stat
%     integer :: mat_type
%     logical :: geo_non_flg
%     integer :: now_thread
%     integer, external :: OMP_GET_THREAD_NUM

%                       try using stack
%     real(8) :: cep_blk_n1(mxvl,nstr,nstr), P_blk_n1(mxvl,nstrs)
%     real(8) :: detF(mxvl),fn1inv(mxvl,ndim,ndim),rnh(mxvl,ndim,ndim)
%     real(8) :: fnh(mxvl,ndim,ndim),dfn(mxvl,ndim,ndim)
%     real(8) :: fnhinv(mxvl,ndim,ndim),ddt(mxvl,nstr)
%     real(8) :: uddt(mxvl,nstr),qnhalf(mxvl,nstr,nstr)
%     real(8) :: qn1(mxvl,nstr,nstr),A_blk_n1(mxvl,nstrs*nstrs)
%     real(8), allocatable, dimension(:,:,:) :: rnh, fnh, dfn, fnhinv
%     real(8), allocatable, dimension(:) :: detF
%     real(8), allocatable, dimension(:,:,:) :: fn1inv, cep_blk_n1
%     real(8), allocatable, dimension(:,:) :: ddt, uddt
%     real(8), allocatable, dimension(:,:,:) :: qnhalf, qn1
%     real(8), allocatable, dimension(:,:) :: P_blk_n1
%     real(8), allocatable, dimension(:,:) :: A_blk_n1
%
%             cauchy stress and rotation matrix @ n+1
%     real(8) :: qtn1(mxvl,nstr,nstr), cs_blk_n1(mxvl,nstr)

%             allocate and initialization
%     allocate( detF(mxvl) )
%     allocate( fn1inv(mxvl,ndim,ndim) )
%     allocate( cep_blk_n1(mxvl,nstr,nstr) )
%     allocate( rnh(mxvl,ndim,ndim), fnh(mxvl,ndim,ndim),
%    &          dfn(mxvl,ndim,ndim), fnhinv(mxvl,ndim,ndim) )
%             warp3d original allocate and initialization
%             nstr = 6; nstrs = 9;
%     allocate( ddt(mxvl,nstr), uddt(mxvl,nstr),
%    &          qnhalf(mxvl,nstr,nstr), qn1(mxvl,nstr,nstr) )
%     allocate( P_blk_n1(mxvl,nstrs) )
%     allocate( A_blk_n1(mxvl,nstrs*nstrs) )

mxvl = 128;
ndim = 3;
nstr = 6;
nstrs = 9;

detF = zeros(mxvl,1);
fn1inv = zeros(mxvl,ndim,ndim);
cep_blk_n1 = zeros(mxvl,nstr,nstr);
rnh = zeros(mxvl,ndim,ndim);
fnh = zeros(mxvl,ndim,ndim);
dfn = zeros(mxvl,ndim,ndim);
fnhinv = zeros(mxvl,ndim,ndim);
ddt = zeros(mxvl,3,3);
uddt = zeros(mxvl,nstr);
qnhalf = zeros(mxvl,nstr,nstr);
qn1 = zeros(mxvl,nstr,nstr);
cs_blk_n1 = zeros(mxvl,3,3);
P_blk_n1 = zeros(mxvl,3,3);
A_blk_n1 = zeros(mxvl,nstrs*nstrs);
%
%     initialize local_work based on elblks(0,blk) and elblks(1,blk)
%
ngp = fftngp;
gpn = 1;
geo_non_flg = 1;
span = elblks(1, blk);
felem = elblks(2, blk);
currmat = matList(felem);
mat_type = matprp(9,currmat);
adaptive_flag = 0;
%
%
%                        initialize local work
%
local_work.dt = tstep;
local_work.blk = blk;
local_work.span = span;
local_work.felem = felem;
local_work.num_int_points = ngp;
local_work.gpn = gpn;
local_work.step = step;
local_work.iter = iter; % global newton iteration (>1)
local_work.geo_non_flg = geo_non_flg;
local_work.is_cohes_elem = 0;
local_work.block_has_nonlocal_solids = 0;
local_work.mat_type = mat_type;
local_work.matnum = currmat;
local_work.segmental = 0;
local_work.number_points = 0;
local_work.curve_set_number = 0;
local_work.fgm_enode_props = 0;
local_work.killed_status_vec = 0;
local_work.iout = out;
local_work.int_order = 1;
local_work.num_enodes = 8;
local_work.num_enode_dof = 3;
local_work.lnelas_vec = 0;
local_work.bbar_flg = 1;
local_work.material_cut_step = 0;
local_work.adaptive_flag = adaptive_flag;
local_work.eps_bbar = 0;
local_work.is_solid_matl = 1;
local_work.is_umat = 0;
local_work.umat_stress_type = 1;
local_work.is_crys_pls = (mat_type == 10);
local_work.temperatures = 297.0D0;
local_work.temperatures_ref = 297.0D0;
%
%     allocate memory for local_work according to material model
%
%     call mm01_set_sizes( info_vector )
%     hist_size = info_vector(1)
%     cep_size = info_vector(2)
%     block_size = span * ngp * cep_size
%     local_work%hist_size_for_blk = hist_size
%     call recstr_allocate( local_work )
local_work.allocate();
local_work.elem_type = 2; % lsdisop element
%
%     grab material parameters, state variables and history
%     due(nodal displacement from n to n+1)
%     ue(nodal displacement from 0 to n)
%     ce_0, cd_n, cd_mid, cd_n1: nodal coordinate
%                 ***hard coded for now***
%
%                hard code material parameters
%     call mm01_hardCoded
setup_mm01_rknstr( span, adaptive, local_work );
%
%                grab global variables to local_work
%
hist_size = local_work.hist_size_for_blk;
for  j = 1:hist_size
    for  i = 1:span
        local_work.elem_hist(i,j) = history_blocks(blk).ptr(j,i);
    end
end
for  j = 1:nstr
    for  i = 1:span
        local_work.ddtse(i,j) = eps_n_blocks(blk).ptr(j,i);
    end
end

for  j = 1:nstr
    for  i = 1:span
        local_work.strain_n(i,j) = eps_n_blocks(blk).ptr(j,i);
    end
end

for  j = 1:nstrs
    for  i = 1:span
        local_work.urcs_blk_n(i,j) = urcs_n_blocks(blk).ptr(j,i);
    end
end
%
%                grab global Fn and Fn1 to local_work
%
for ii = 1:span
    currElem = felem + ii - 1;
    local_work.fn(ii,1,1)= Fn(currElem,1);
    local_work.fn(ii,1,2)= Fn(currElem,2);
    local_work.fn(ii,1,3)= Fn(currElem,3);
    local_work.fn(ii,2,1)= Fn(currElem,4);
    local_work.fn(ii,2,2)= Fn(currElem,5);
    local_work.fn(ii,2,3)= Fn(currElem,6);
    local_work.fn(ii,3,1)= Fn(currElem,7);
    local_work.fn(ii,3,2)= Fn(currElem,8);
    local_work.fn(ii,3,3)= Fn(currElem,9);
end
for ii = 1:span
    currElem = felem + ii - 1;
    local_work.fn1(ii,1,1)= Fn1(currElem,1);
    local_work.fn1(ii,1,2)= Fn1(currElem,2);
    local_work.fn1(ii,1,3)= Fn1(currElem,3);
    local_work.fn1(ii,2,1)= Fn1(currElem,4);
    local_work.fn1(ii,2,2)= Fn1(currElem,5);
    local_work.fn1(ii,2,3)= Fn1(currElem,6);
    local_work.fn1(ii,3,1)= Fn1(currElem,7);
    local_work.fn1(ii,3,2)= Fn1(currElem,8);
    local_work.fn1(ii,3,3)= Fn1(currElem,9);
end
for ii = 1:span
    fnh(ii,1,1) = half * ( local_work.fn(ii,1,1) ...
        &              + local_work.fn1(ii,1,1) );
    fnh(ii,1,2) = half * ( local_work.fn(ii,1,2) ...
        &              + local_work.fn1(ii,1,2) );
    fnh(ii,1,3) = half * ( local_work.fn(ii,1,3) ...
        &              + local_work.fn1(ii,1,3) );
    fnh(ii,2,1) = half * ( local_work.fn(ii,2,1) ...
        &              + local_work.fn1(ii,2,1) );
    fnh(ii,2,2) = half * ( local_work.fn(ii,2,2) ...
        &              + local_work.fn1(ii,2,2) );
    fnh(ii,2,3) = half * ( local_work.fn(ii,2,3) ...
        &              + local_work.fn1(ii,2,3) );
    fnh(ii,3,1) = half * ( local_work.fn(ii,3,1) ...
        &              + local_work.fn1(ii,3,1) );
    fnh(ii,3,2) = half * ( local_work.fn(ii,3,2) ...
        &              + local_work.fn1(ii,3,2) );
    fnh(ii,3,3) = half * ( local_work.fn(ii,3,3) ...
        &              + local_work.fn1(ii,3,3) );
end
dfn = local_work.fn1 - local_work.fn;
%
%              compute rotation tensor
%
for ii = 1:span
    rnh(ii,:,:) = poldec(fnh(ii,:,:));
    rn1 = poldec(local_work.fn1(ii,:,:));
    local_work.rot_blk_n1(ii,:) = reshape(rn1,1,9);
end
%
%              compute uddt( detF is a dummy argument )
%
for ii = 1:span
    ddt(ii,:,:) = dfn(ii,:,:)/fnh(ii,:,:);
    uddt_temp = rnh(ii,:,:)' * ddt(ii,:,:) * rnh(ii,:,:);
    uddt(ii,1) = uddt_temp(1,1);
    uddt(ii,2) = uddt_temp(2,2);
    uddt(ii,3) = uddt_temp(3,3);
    uddt(ii,4) = uddt_temp(1,2) + uddt_temp(2,1);
    uddt(ii,5) = uddt_temp(2,3) + uddt_temp(3,2);
    uddt(ii,6) = uddt_temp(1,3) + uddt_temp(3,1);
end
%
%     recover stress and stiffness
%
rstgp1( local_work, uddt, cep_blk_n1 );
%
%           For CP model, calculate the gradient of the elastic
%           rotations at the element level by linear curve fit.
%
%           For linear models this will just be based on the plastic rotations
%           which may or may not be a realistic assumption
%
%     if( local_work.mat_type == 10 )
%         call rknstr_finish_cp
%     end
%
%     scatter local variables to global
%
rplstr( span, felem, ngp, mat_type, iter, geo_non_flg, local_work, blk, ...
urcs_n1_blocks, eps_n1_blocks, rot_n1_blocks, history1_blocks);
%
%     compute cs_blk_n1(Cauchy stress) from urcs_blk_n1
%
temp = [1,4,6;4,2,5;6,5,3];
for ii = 1:span
    rn1 = reshape(local_work.rot_blk_n1(ii,:),3,3);
    urcsn1_v = local_work.urcs_blk_n1(ii,1:6);
    urcsn1 = urcsn1_v(temp);
    cs_blk_n1(ii,1:3,1:3) = rn1 * urcsn1 * rn1';
end
%
%     pull back cs_blk_n1(Cauchy stress) to 1st PK stress
%     P = J*sigma*F^{-T}
for ii = 1:span
    f_local = local_work.fn1(ii,1:3,1:3);
    cs_local = cs_blk_n1(ii,1:3,1:3);
    p_local = det(f_local) * cs_local / (f_local');
    P_blk_n1(ii,1:9) = reshape(p_local',1,9);
end
%
%         rotate from ( d_urcs / d_uddt ) to ( Green-Naghdi / D )
%     (1) extract stiffness from mod_eleblocks
%     (2) push forward to current configuration
%     (3) store in cep_blk_n1(mxvl,nstr,nstr)
%
for ii = 1:span
    % take cep from cep_blk_n1
    temp = (ii-1)*21;
    cepn1 = zeros(6,6);
    cepn1(1,1) = cep_blk_n1(blk).ptr(temp+1);
    cepn1(2,1:2) = cep_blk_n1(blk).ptr(temp+2:temp+3);
    cepn1(3,1:3) = cep_blk_n1(blk).ptr(temp+4:temp+6);
    cepn1(4,1:4) = cep_blk_n1(blk).ptr(temp+7:temp+10);
    cepn1(5,1:5) = cep_blk_n1(blk).ptr(temp+11:temp+15);
    cepn1(6,1:6) = cep_blk_n1(blk).ptr(temp+16:temp+21);
    
    % form 6x6 matrix
    % push forward to current configuration
    % pull back to reference configuration
end
gptns1( local_work, cep_blk_n1, qtn1 );
%
%     pull back cep_blocks to dP/dF
%     assume that Green-Naghdi rate is close to Lie derivative
%     see Simo & Hughes, chapter 7
%     input:  cep_blk_n1(mxvl,nstr,nstr)
%     output: A_blk_n1(mxvl,nstrs*nstrs)
call cep2A( span, cs_blk_n1, cep_blk_n1, fn1inv, detF, A_blk_n1, out);
%
%     update global P and tangent stiffness
for ii = 1:span
    currElem = felem + ii - 1;
    Pn1(currElem, 1:9) = P_blk_n1(ii, 1:9);
    K4(currElem, 1:81) = A_blk_n1(ii, 1:81);
end
%
%     get cut-step flag and update plastic work
%
%     do ii = 1, span
%       currElem = felem + ii - 1
%       call constitutive(Fn1(currElem,:), matList(currElem),
%    &                     Pn1(currElem,:), K4(currElem,:))
%     enddo
%
% 100 deallocate( rnh, fnh, dfn, fnhinv, fn1inv )
%     deallocate( ddt, uddt, qnhalf, qn1 )
%     deallocate( P_blk_n1 )
%     deallocate( A_blk_n1 )
end
%
%
%     ****************************************************************
%     *                                                              *
%     *                   subroutine setup_mm01_rknstr               *
%     *                                                              *
%     *                       written by : RM                        *
%     *                                                              *
%     *                   last modified : 5/17/2018 RM               *
%     *                                                              *
%     *     set up material model #1 (bilinear mises) for stress     *
%     *     updating                                                 *
%     *                                                              *
%     ****************************************************************
%
function setup_mm01_rknstr( span, adaptive, local_work )
%     use fft, only: matprp, lmtprp, imatprp, dmatprp, smatprp
%     implicit none
%     include 'param_def'
%     include 'include_sig_up'

%                     global
%     integer :: span
%     logical :: adaptive

%                     local
%     integer :: currmat, i
%     real(8) :: ym, nu, beta, tan_e, yld, alpha, rho

currmat = local_work.matnum;

for i = 1:span
    ym                       = matprp(1,currmat);
    nu                       = matprp(2,currmat);
    yld                      = matprp(5,currmat);
    tan_e                    = matprp(4,currmat);
    rho                      = matprp(7,currmat);
    alpha                    = matprp(6,currmat);
    beta                     = matprp(3,currmat);
    local_work.e_vec(i)      = ym;
    local_work.nu_vec(i)     = nu;
    local_work.beta_vec(i)   = beta;
    local_work.tan_e_vec(i)  =  tan_e;
    local_work.sigyld_vec(i) = yld;
    local_work.h_vec(i)      = tan_e*ym/(ym - tan_e) ;
    local_work.e_vec_n(i)    = ym;
    local_work.nu_vec_n(i)   = nu;
end

end