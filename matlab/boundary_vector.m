function [B,Kb] = boundary_vector_bis(p,e,variables)
%BOUNDARY_VECTOR Outputs boundary vector
%   Henri De Plaen, KU Leuven
%
%  Input:
%         p : Node coordinates
%         e : Edge vector (bound1: index 1-5, bound2: index 6)
%         t : Triangle vertices
%
%  Output:
%       B  : Boundary vector
%            sparse np

%% FUNC BOUNDARY
np = size(p,1) ;

e6 = e(e(:,5)==6,:) ;
e1 = e6(:,1) ; 
e2 = e6(:,2) ;

hp = p(e2,:) - p(e1,:) ; 
h = sqrt(sum(hp.^2,2)) ;

[bu,bv,kbu,kbv] = buv(p,variables) ;

bbu = h.*(bu(e1)+bu(e2))/2 ;
bbv = h.*(bv(e1)+bv(e2))/2 ;
kbbu = h.*(kbu(e1)+kbu(e2))/2 ;
kbbv = h.*(kbv(e1)+kbv(e2))/2 ;

B = sparse(e1,1,bbu,2*np,1)+sparse(e2,1,bbu,2*np,1) + ...
    sparse(e1+np,1,bbv,2*np,1)+sparse(e2+np,1,bbv,2*np,1) ;

Kb = sparse(e1,e1,kbbu,2*np,2*np)+sparse(e2,e2,kbbu,2*np,2*np) + ...
    sparse(e1+np,e1+np,kbbv,2*np,2*np)+sparse(e2+np,e2+np,kbbv,2*np,2*np) ;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bu,bv,kbu,kbv] = buv(p,variables)
%% BOUNDARY VECTOR
r = p(:,1) ;
np = size(p,1) ;

% bu = r*-variables.hu*variables.Cuamb ;
% bv = r*variables.hv*variables.Cvamb ;
% 
% kbu = variables.hu*r ;
% kbv = variables.hv*r ;

bu = -variables.hu*variables.Cuamb.*ones(np,1) ;
bv = -variables.hv*variables.Cvamb.*ones(np,1) ;

kbu = variables.hu.*ones(np,1) ;
kbv = variables.hv.*ones(np,1) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%