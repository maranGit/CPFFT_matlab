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
function [nelblk, elblks] = inelbk( matList, noelem )
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
      param_def
      auto_size    = mxvl;  %  from param_def
      elblks = zeros(4,mxnmbl);
%     display      = 1;
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
      elblks(2,1)  = 1;  % first element in block
      elblks(1,1)  = 1;  % number elements in block
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
         if( current_size > auto_size )
             newblk = 1;
         end
         if( ~ newblk )
           ele_matmodel = matList( element );
           compatible = ( blk_matmodel == ele_matmodel );
           if( ~ compatible ) 
               newblk = 1;
           end
         end
         if( ~ newblk )
           elblks(1,nelblk)  = current_size;
           continue
         end
         nelblk = nelblk + 1;
         if( nelblk > mxnmbl )
%           param = nelblk;
%      call errmsg(74,param,dums,dumr,dumd)
%      call die_abort
            error('too many element blocks required');
         end
         felem             = element;
         elblks(2,nelblk)  = felem;
         elblks(1,nelblk)  = 1;
         current_size      = 1;
         blk_matmodel = matList( felem );
      end % on element

end