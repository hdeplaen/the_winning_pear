function [Bu,Bv,Kbu,Kbv] = boundary_vector(p,e,variables,refer)
%Outputs boundary vector
%   [B,Kb] = boundary_vector(p,e,variables,axis_e,refer)
%
%   Henri De Plaen, KU Leuven
%
%  Input:
%         p : Node coordinates
%         e : Edge vector (bound1: index 1-5, bound2: index 6)
%         t : Triangle vertices
%         refer : 'euclid' or 'radial' (legacy)
%
%  Output:
%       Bu  : Independent boundary vector term
%             sparse of Cu (np,1)
%       Bv  : Independent boundary vector term
%             sparse of Cv (np,1)
%       Kbu : Boundary vector term dependent on Cu
%            sparse (np,np)
%       Kbv : Boundary vector term dependent on Cv
%            sparse (np,np)

%% FUNC BOUNDARY
np = size(p,1) ;                                % number of points (boundary + non-boundary)

g1 = e(:,3) == 1 ;                              % index of gamma1 boundary segments (on the curve, not the symmetry axis)

e1 = e(g1,1) ;                                  % index of first points of gamma1 segments
e2 = e(g1,2) ;                                  % index of second points of gamma1 segements

hp = p(e2,:) - p(e1,:) ;
h = sqrt(sum(hp.^2,2)) ;                        % length of the segments

[bu,bv,kbu,kbv] = buv(p,variables, refer) ;     % estimation of the boundary value at the nodes

bbu1 = h.*(2*bu(e1)+bu(e2))/6 ;                 % estimation of the boundary vector values for each segment
bbu2 = h.*(bu(e1)+2*bu(e2))/6 ;

bbv1 = h.*(2*bv(e1)+bv(e2))/6 ;
bbv2 = h.*(bv(e1)+2*bv(e2))/6 ;

kbbud1 = h.*(3*kbu(e1)+kbu(e2))/12 ;
kbbud2 = h.*(kbu(e1)+3*kbu(e2))/12 ;
kbbum  = h.*(kbu(e1)+kbu(e2))/12 ;

kbbvd1 = h.*(3*kbv(e1)+kbv(e2))/12 ;
kbbvd2 = h.*(kbv(e1)+3*kbv(e2))/12 ;
kbbvm  = h.*(kbv(e1)+kbv(e2))/12 ;

% construction of the independent boundary vector
Bu = sparse(e1,1,bbu1,np,1) + sparse(e2,1,bbu2,np,1) ;
Bv = sparse(e1,1,bbv1,np,1) + sparse(e2,1,bbv2,np,1) ;

% construction of the linearly dependant (on Cu and Cv) boundary vector
Kbu = sparse(e1,e1,kbbud1,np,np) + sparse(e2,e2,kbbud2,np,np) + ...
      sparse(e1,e2,kbbum,np,np) + sparse(e2,e1,kbbum,np,np) ;
Kbv = sparse(e1,e1,kbbvd1,np,np) + sparse(e2,e2,kbbvd2,np,np) + ...
      sparse(e1,e2,kbbvm,np,np) + sparse(e2,e1,kbbvm,np,np) ;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bu,bv,kbu,kbv] = buv(p,variables, refer)
% Estimation of the boundary value at the given nodes

%% BOUNDARY VECTOR
r = p(:,1) ;

switch refer
    case 'radial'
        % independent boundary terms
        bu = r*variables.hu*variables.Cuamb ;
        bv = r*variables.hv*variables.Cvamb ;
        
        % boundary terms dependent on Cu and Cv
        kbu = -variables.hu*r ;
        kbv = -variables.hv*r ;
    case 'cartesian'
        np = size(p,1) ;  % number of points
        
        % independent boundary terms
        bu = variables.hu*variables.Cuamb.*ones(np,1) ;
        bv = variables.hv*variables.Cvamb.*ones(np,1) ;
        
        % boundary terms dependent on Cu and Cv
        kbu = -variables.hu.*ones(np,1) ;
        kbv = -variables.hv.*ones(np,1) ;
    otherwise ; error('Wrong referential') ;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%