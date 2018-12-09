function int_function_matrix = int_function(p,t)
%Outputs function vector for linearization F(x) + dF(x)(x1-x0)
%   [F,dF] = fun_vector(p,t,sol,variables)
%
%  Input:
%         p : Node coordinates
%         t : Triangle coordinates
%
%  Output:
%       int_function_matrix : Sparse matrix used to interpolate source
%                             term (function Ru/Rv)
%
% Henri De Plaen, KU Leuven

np = size(p,1) ;
nt = size(t,1) ;

int_function_matrix = zeros(np,np) ; % prealloc

for idx = 1:nt
    loc_t   = t(idx,1:3) ;
    r       = p(loc_t,1) ;
    z       = p(loc_t,2) ;
    int_function_matrix(loc_t,loc_t) = int_function_matrix(loc_t,loc_t) +...
        int_function_block(r,z);
end

int_function_matrix = sparse(int_function_matrix) ;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function int_function_matrix_block = int_function_block(r,z)

T = triangle_areas(r,z);       
       int_function_matrix_block = T./60.* ...
           [[6*r(1)+2*r(2)+2*r(3) 2*r(1)+2*r(2)+r(3) 2*r(1)+r(2)+2*r(3)] ;
           [2*r(1)+2*r(2)+r(3) 2*r(1)+6*r(2)+2*r(3) r(1)+2*r(2)+2*r(3)] ;
           [2*r(1)+r(2)+2*r(3) r(1)+2*r(2)+2*r(3) 2*r(1)+2*r(2)+6*r(3)]] ;
       
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T] = triangle_areas(r,z)

    T = r(2)*z(3)+r(1)*z(2) + r(3)*z(1) - r(2)*z(1) - r(1)*z(3) - r(3)*z(2) ; 
    
end