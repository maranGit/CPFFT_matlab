classdef blocks_ptr_type < handle % pass by reference
    properties
        ptr
    end
    
    properties(Constant)
        mxvl = 128;
        nstr = 6;
        nstrs = 9;
    end
    
    methods
        function obj = blocks_ptr_type()
            
        end
        
        function allocate(obj,dim1)
            obj.ptr = zeros(dim1,1);
        end
    end
end