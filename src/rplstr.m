%     ****************************************************************
%     *                                                              *
%     *                      subroutine rplstr                       *
%     *                                                              *
%     *                       written by : rhd                       *
%     *                                                              *
%     *                   last modified : 7/28/2016 rhd              *
%     *                                                              *
%     *     stores globally the recovered material                   *
%     *     and stress states.                                       *
%     *                                                              *
%     ****************************************************************
function rplstr( span, felem, ngp, mat_type, iter, geonl, local_work, blk, ...
    urcs_n1_blocks, eps_n1_blocks, rot_n1_blocks, history1_blocks)
%
%
%             replace stress/strain data:
%
%               1) unrotated cauchy stresses at n+1 (all models)
%
%               2) ddtse comes back from element strain routine with
%                  total strains at n+1 (all models)
%
%               3) [R,n+1] for geonl models
%
%               4) history data updated to n+1
%
%             iter = 0 indicates the "pre" step computations used
%             by the load step processor to get estimated stresses
%             for non-zero imposed displacements and/or the
%             extrapolated displacement increment for step.
%             we do no want to update the material state variables
%             especially (2) above since they most likely will be
%             used immediately after this for stiffness
%             computation. Except for the CP model since it hides the
%             [D] matrix in there.
%
%
process_hist = mat_type == 1  ||  mat_type == 2 || ...
    mat_type == 3  ||  mat_type == 4 || ...
    mat_type == 5  ||  mat_type == 6 || ...
    mat_type == 7  ||  mat_type == 8 || ...
    mat_type == 10;
%
for j = 1:nstrs
    for i = 1:span
        temp = (j-1)*nstrs + i;
        urcs_n1_blocks(blk).ptr(temp) = local_work.urcs_blk_n1(i,j);
    end
end
%
for j = 1:nstr
    for i = 1:span
        temp = (j-1)*nstr + i;
        eps_n1_blocks(blk).ptr(temp) = local_work.ddtse(i,j);
    end
end
%
if( iter > 0 && geonl )
    for j = 1:9
        for i = 1:span
            temp = (j-1)*nstrs + i;
            rot_n1_blocks(blk).ptr(temp) = local_work.rot_blk_n1(i,j);
        end
    end
end
%
save_history_1 = process_hist && iter > 0;
save_history_2 = mat_type == 10; % [D] is hidden in history ...
if( save_history_1  ||  save_history_2 )
    hist_size = history_blk_list(blk);
    for j = 1:hist_size
        for i = 1:span
            temp = (j-1)*hist_size + i;
            history1_blocks(blk).ptr(temp) = local_work.elem_hist1(i,j);
        end
    end
end
end