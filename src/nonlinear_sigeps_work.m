classdef nonlinear_sigeps_work < handle
    properties
        ce_0;                % (mxvl,mxecor)
        ce_n;                % (mxvl,mxecor)
        ce_mid;              % (mxvl,mxecor)
        ce_n1;               % (mxvl,mxecor)
        trnmte;              % (mxvl,mxedof,mxndof)
        det_j;               % (mxvl,mxgp)
        det_j_mid;           % (mxvl,mxgp)
        weights;             % mxgp
        nxi;                 % (mxndel,mxgp)
        neta;                % (mxndel,mxgp)
        nzeta;               % (mxndel,mxgp)
        gama;                % (mxvl,3,3,mxgp)
        gama_mid;            % (mxvl,3,3,mxgp)
        fn;                  % (mxvl,3,3)
        fn1;                 % (mxvl,3,3)
        dfn1;                % (mxvl)
        vol_block;           % (mxvl,8,3) )
        volume_block;        % (mxvl)
        volume_block_0;      % (mxvl)
        volume_block_n;      % (mxvl)
        volume_block_n1;     % (mxvl)
        jac;                 % (mxvl,3,3)
        b;                   % (mxvl,mxedof,nstr)
        ue;                  % (mxvl,mxedof)
        due;                 % (mxvl,mxedof)
        uenh;                % (mxvl,mxedof)
        uen1;                % (mxvl,mxedof)
        urcs_blk_n;          % (mxvl,nstrs,mxgp)
        urcs_blk_n1;         % (mxvl,nstrs,mxgp)
        rot_blk_n1;          % (mxvl,9,mxgp)
        rtse;                % (mxvl,nstr,mxgp)
        elem_hist1;
        elem_hist;
        ddtse;               % (mxvl,nstr,mxgp)
        strain_n;            % (mxvl,nstr,mxgp)
        dtemps_node_blk;     % (mxvl,mxndel)
        temps_ref_node_blk;  % (mxvl,mxndel)
        temps_node_blk;      % (mxvl,mxndel)
        temps_node_ref_blk;  % (mxvl,mxndel)
        cohes_temp_ref;      % (mxvl)
        cohes_dtemp;         % (mxvl)
        cohes_temp_n;        % (mxvl)
        nu_vec;              % (mxvl)
        beta_vec;            % (mxvl)
        h_vec;               % (mxvl)
        e_vec;               % (mxvl)
        sigyld_vec;          % (mxvl)
        alpha_vec;           % (mxvl,6)
        e_vec_n;             % (mxvl)
        nu_vec_n;            % (mxvl)
        gp_sig_0_vec;        % (mxvl)
        gp_sig_0_vec_n;      % (mxvl)
        gp_h_u_vec;          % (mxvl)
        gp_h_u_vec_n;        % (mxvl)
        gp_beta_u_vec;       % (mxvl)
        gp_beta_u_vec_n;     % (mxvl)
        gp_delta_u_vec;      % (mxvl)
        gp_delta_u_vec_n;    % (mxvl)
        alpha_vec_n;         % (mxvl,6)
        h_vec_n;             % (mxvl)
        n_power_vec;         % (mxvl)
        f0_vec;              % (mxvl)
        eps_ref_vec;         % (mxvl)
        m_power_vec;         % (mxvl)
        q1_vec;              % (mxvl)
        q2_vec;              % (mxvl)
        q3_vec;              % (mxvl)
        nuc_s_n_vec;         % (mxvl)
        nuc_e_n_vec;         % (mxvl)
        nuc_f_n_vec;         % (mxvl)
        
        %
        eps_curve;           % (max_seg_points)
        shape;               % (mxndel,mxgp)
        characteristic_length; % (mxvl)
        intf_prp_block;      % (mxvl,50)
        cohes_rot_block;     % (mxvl,3,3)
        enode_mat_props;     % (mxndel,mxvl,mxndpr)
        tan_e_vec;           % (mxvl)
        fgm_flags;           % (mxvl,mxndpr)
        mm05_props;          % (mxvl,10)
        mm06_props;          % (mxvl,10)
        mm07_props;          % (mxvl,10)
        umat_props;          % (mxvl,50)
        
        dt;
        total_model_time;
        time_n;
        beta_fact;
        block_energy;
        eps_bbar;
        block_plastic_work;
        step_scale_fact;
        alpha_dmg;
        ls;
        ll;
        lt;
        sv;
        lv;
        tv;      % 3 each
        %
        felem;
        elem_type;
        matnum;
        int_order;
        mat_type;
        num_enodes;
        num_enode_dof;
        totdof;
        num_int_points;
        span;
        iter;
        step;
        gpn;
        number_points;
        cohes_type;
        curve_set_number;
        surface;
        hist_size_for_blk;
        iout;
        blk;
        umat_stress_type;
        cep_sym_size;
        num_threads;
        inter_mat;
        macro_sz;
        cp_sz;
        now_thread;
        %
        trne;                     % (mxvl,mxndel)
        geo_non_flg;
        bbar_flg;
        trn_e_block;
        trn_e_flags;
        first;
        material_cut_step;
        signal_flag;
        adaptive_flag;
        temperatures;
        lnelas_vec;
        nuc_vec;
        nonlinear_flag;
        allow_cut;
        segmental;
        power_law;
        temps_node_to_process;
        temperatures_ref;
        fgm_enode_props;
        is_cohes_elem;
        linear_displ_elem;
        adjust_const_elem;
        is_axisymm_elem;
        killed_status_vec;
        block_killed;
        is_umat;
        is_solid_matl;
        is_crys_pls;
        compute_f_bar;
        compute_f_n;
        is_cohes_nonlocal;
        is_inter_dmg;
        block_has_nonlocal_solids;
        %
        %                       nonlocal cohesive support
        %
        %    &            top_surf_solid_stresses_n(:,:),             % (mxvl,nstrs)
        %    &            bott_surf_solid_stresses_n(:,:),            % (mxvl,nstrs)
        %    &            top_surf_solid_eps_n(:,:),                  % (mxvl,nstr)
        %    &            bott_surf_solid_eps_n(:,:),                 % (mxvl,nstr)
        %    &            nonlocal_stvals_bott_n(:,:),                % (mxvl,nxx)
        %    &            nonlocal_stvals_top_n(:,:)                  % (mxvl,nxx)
        %
        %              workspace to store block of solid matl. nonlocal values
        %              while processing an integration point
        %
        % % nonlocal cohesive support
        %    &            nonlocal_state_blk(:,:)                     % (mxvl,nxx)
        %
        %             where nxx above = nonlocal_shared_state_size
        %
        %      integer, allocatable :: top_surf_solid_elements(:),    % (mxvl)
        %    &                         bott_surf_solid_elements(:)
        %      integer, allocatable :: top_solid_matl(:),             % (mxvl)
        %    &                         bott_solid_matl(:),
        %    &                         nstacks(:), nper(:)
        %
        %                      Added stuff for CP
        %
        %     logical, allocatable :: debug_flag(:)                   % mxvl
        %     double precision, allocatable :: local_tol(:)           % mxvl
        %     integer, allocatable :: ncrystals(:)                    % mxvl
        %     integer, allocatable :: angle_type(:)                   % mxvl
        %     integer, allocatable :: angle_convention(:)             % mxvl
        %     type(crystal_properties), allocatable :: c_props(:,:)   % mxvl,max_crystals
    end
    
    methods
        function obj = nonlinear_sigeps_work()
        end
        function allocate(local_work)
            mxvl = 128;
            nstrs = 9;
            nstr = 6;
            mxndel = 20;
            
            local_work.fn = zeros(mxvl,3,3);
            local_work.fn1 = zeros(mxvl,3,3);
            local_work.dfn1 = zeros(mxvl);
            
            local_work.urcs_blk_n = zeros(mxvl,nstrs);
            local_work.urcs_blk_n1 = zeros(mxvl,nstrs);
            local_work.rot_blk_n1 = zeros(mxvl,9);
            local_work.rtse = zeros(mxvl,nstr);
            %
            local_work.ddtse = zeros(mxvl,nstr);
            local_work.strain_n = zeros(mxvl,nstr);
            local_work.dtemps_node_blk = zeros(mxvl,mxndel);
            local_work.temps_ref_node_blk = zeros(mxvl,mxndel);
            local_work.temps_node_blk = zeros(mxvl,mxndel);
            local_work.temps_node_ref_blk = zeros(mxvl,mxndel);
            local_work.nu_vec = zeros(mxvl);
            local_work.beta_vec = zeros(mxvl);
            local_work.h_vec = zeros(mxvl);
            local_work.tan_e_vec = zeros(mxvl);
            local_work.e_vec = zeros(mxvl);
            %
            local_work.sigyld_vec = zeros(mxvl);
            local_work.alpha_vec = zeros(mxvl,6);
            local_work.e_vec_n = zeros(mxvl);
            local_work.nu_vec_n = zeros(mxvl);
            %
            local_work.alpha_vec_n = zeros(mxvl,6);
            local_work.h_vec_n = zeros(mxvl);
            local_work.n_power_vec = zeros(mxvl);
            %
            local_mt = local_work.mat_type;
            %
            span_local                   = local_work.span;
            blk_local                    = local_work.blk;
            ngp                          = local_work.num_int_points;
            hist_size                    = history_blk_list(blk_local);
            local_work.hist_size_for_blk = hist_size;
            %
            local_work.elem_hist1 = zeros(span_local,hist_size);
            local_work.elem_hist = zeros(span_local,hist_size);
            %
            local_work.weights = 0;
        end
    end
end