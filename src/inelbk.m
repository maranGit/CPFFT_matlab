%
%     ****************************************************************
%     *                                                              *
%     *                 subroutine      inelbk                       *
%     *                                                              *
%     *                       written by : RM                        *
%     *                                                              *
%     *                   last modified: 1/22/18                     *
%     *                                                              *
%     *     automatic generation of element blocking: simple version *
%     *     for threaded w/o vectorized blocking. optional           *
%     *     assignment of blocks to domains w/ simple algorithm      *
%     *                                                              *
%     ****************************************************************
%
function inelbk( matList )
%     implicit none
%     include 'common.main'
%
%     integer :: auto_size
%     logical :: display
%     integer :: felem, current_size, element, param, i
%     logical :: newblk, compatible
%     integer :: blk_matmodel, ele_matmodel
%     integer, intent(in) :: matList(noelem)
%
%     Ran: hard coded variables for autoblock
%
      auto_size    = mxvl;  %  from param_def
      display      = 1;
%
%                     first generation of automatic assignment of
%                     elements to blocks.
%
%                     a) sequential pass thru elements
%                     b) assign to current block or open new
%                        block if full or element/material
%                        combinations are not compatible with
%                        current block.
%                     c) no element renumbering or vectorized
%                        blocking in this first set of
%                        features.
%
%
      nelblk       = 1;   % in common main
      current_size = 1;
      felem        = 1;
      elblks(1,1)  = 1;  % first element in block
      elblks(0,1)  = 1;  % number elements in block
      if( noelem == 1 ) 
          return;
      end
%
      blk_matmodel = matList( 1 );
%
%                     1. element fits in current block?
%                     2. if yes, is it compatible with elements now
%                        in the block?
%                     3. if yes, update current number of elements
%                        in the block, next element
%                     4. otherwise start a new block. set first
%                        element in the block, init block size,
%                        load props for first element in block for
%                        subsequent comparisons
%
      for element = 2:noelem
         newblk       = 0;
         current_size = current_size + 1;
         if( current_size .gt. auto_size )
             newblk = 1;
         end
         if( ~ newblk )
           ele_matmodel = matList( element );
           compatible = ( blk_matmodel .eq. ele_matmodel );
           if( ~ compatible ) 
               newblk = 1;
           end
         end
         if( ~ newblk )
           elblks(0,nelblk)  = current_size;
           cycle;
         end
         nelblk = nelblk + 1;
         if( nelblk .gt. mxnmbl )
            param = nelblk;
%      call errmsg(74,param,dums,dumr,dumd)
%      call die_abort
            error('too many element blocks required');
         end
         felem             = element;
         elblks(1,nelblk)  = felem;
         elblks(0,nelblk)  = 1;
         current_size      = 1;
         blk_matmodel = matList( felem );
      end % on element

end