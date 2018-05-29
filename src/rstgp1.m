%     ****************************************************************
%     *                                                              *
%     *                      subroutine rstgp1                       *
%     *                                                              *
%     *                       written by : rhd                       *
%     *                                                              *
%     *                   last modified : 1/2/2016 rhd               *
%     *                                                              *
%     *     supervise the computation of strains, stresses and       *
%     *     accompaning stress data at an integration point          *
%     *     for a block of similar elements that use the same        *
%     *     material model code                                      *
%     *                                                              *
%     *              ** geometric nonlinear version **               *
%     *                                                              *
%     ****************************************************************
%
%
%      subroutine rstgp1( props, lprops, iprops, local_work )
function rstgp1( local_work, uddt, cep_blk_n1 )
%
%        allocate and zero. only span rows are used but
%        array operators (e.g uddt = ..) will operate on
%        full content and access uninitialized values.
%
%     allocate( ddt(mxvl,nstr), uddt(mxvl,nstr),
%    &          qnhalf(mxvl,nstr,nstr), qn1(mxvl,nstr,nstr) )

drive_01_update( local_work, uddt, cep_blk_n1 );
end

function drive_01_update( local_work, uddt, cep_blk_n1 )
%      use segmental_curves, only : max_seg_points
%     use elem_block_data, only : gbl_cep_blocks => cep_blocks
%      use main_data, only : extrapolated_du, non_zero_imposed_du
%
%     implicit none
%     include 'param_def'
%     integer, parameter :: max_seg_points=20
%
%                      parameter declarations
%
%      real ::    props(mxelpr,*)   ! all same but read only
%      logical :: lprops(mxelpr,*)
%      integer :: iprops(mxelpr,*)
%     integer :: gpn, iout
%     double precision ::  uddt_displ(mxvl,nstr)
%     include 'include_sig_up'
%
%                       locally defined variables
%
%     integer :: span, felem, type, order, ngp, nnode, ndof, step,
%    &           iter, now_blk, mat_type, number_points, curve_set,
%    &           hist_size_for_blk, curve_type, elem_type, i
%
%     double precision ::
%    &  dtime, gp_temps(mxvl), gp_rtemps(mxvl), gp_dtemps(mxvl),
%    &  zero,  ddummy(1), gp_alpha, ymfgm, et, uddt_temps(mxvl,nstr),
%    &  uddt(mxvl,nstr), cep(mxvl,6,6)
%
%     logical :: geonl, local_debug, temperatures, segmental,
%    &           temperatures_ref, fgm_enode_props
%
%     data zero / 0.0d0 /
%
%           vectorized mises plasticity model with constant hardening
%           modulus. the model supports temperature dependence of
%           the elastic modulus, nu, hprime, and thermal
%           expansion alpha can vary. temperature dependent
%           properties enter through segmental curves.
%
span              = local_work.span;
now_blk           = local_work.blk;
%
%            now standard update process. use nonlinear update and [D]
%            computation for iter = 0 and extrapolation or iter > 1
%            for iter = 0 and no extrapolation, use linear-elastic [D]
%            with props at n+1.
%      if( iter >= 1 .or. extrapolated_du ) then !nonlinear update
%
%     ****************************************************************
%     *                                                              *
%     *                 Ran modify this line                         *
%     *                 no need to call drive_01_update_a            *
%     *                                                              *
%     ****************************************************************
%     if( iter >= 1 ) then !nonlinear update
for ii = 1:span
    ym_n1 = local_work.e_vec(ii);
    nu_n1 = local_work.nu_vec(ii);
    beta = local_work.beta_vec(ii);
    hprime_n1 = local_work.h_vec(ii);
    yld_n1 = local_work.sigyld_vec(ii);
    cgn = local_work.urcs_blk_n(ii,1:6);
    uddt_local = uddt(ii,1:6);
    currhist = local_work.elem_hist(ii,1:11);
    ym_n = ym_n1;
    nu_n = nu_n1;
    [cgn1,currhist1,rtse,yield] = ...
        mm01(ym_n1,nu_n1,beta,hprime_n1,yld_n1,cgn,uddt_local,currhist,ym_n,nu_n);
    cep = cnst1(rtse,nu_n1,ym_n1,currhist1(2),currhist1(5),beta,currhist1(1),1,1,yield);
    local_work.urcs_blk_n1(ii,1:6) = cgn1;
    local_work.elem_hist1(ii,1:11) = currhist1;
    temp = (ii-1) * 21;
    cep_blk_n1(now_blk).ptr(temp+1) = cep(1,1);
    cep_blk_n1(now_blk).ptr(temp+2:temp+3) = cep(2,1:2);
    cep_blk_n1(now_blk).ptr(temp+4:temp+6) = cep(3,1:3);
    cep_blk_n1(now_blk).ptr(temp+7:temp+10) = cep(4,1:4);
    cep_blk_n1(now_blk).ptr(temp+11:temp+15) = cep(5,1:5);
    cep_blk_n1(now_blk).ptr(temp+16:temp+21) = cep(6,1:6);
end
%
%
end