function [B,Kb] = boundary_vector(p,e,variables,axis_e,refer)
%Outputs boundary vector
%   [B,Kb] = boundary_vector(p,e,variables,axis_e,refer)
%
%   Henri De Plaen, KU Leuven
%
%  Input:
%         p : Node coordinates
%         e : Edge vector (bound1: index 1-5, bound2: index 6)
%         t : Triangle vertices
%         axis_e : part of the boundary that corresponds to the axis
%         refer: 'euclid' or 'radial'
%
%  Output:
%       B  : Independent boundary vector term
%            sparse (2*np,1)
%       Kb : Boundary vector term dependent on Cu and Cv
%            sparse (2*np,2*np)

%% FUNC BOUNDARY
np = size(p,1) ;                        % number of points (boundary + non-boundary)

g1 = e(:,5) ~= axis_e ;                 % index of gamma1 boundary segments (on the curve, not the symmetry axis)

e1 = e(g1,1) ;                          % index of first points of gamma1 segments
e2 = e(g1,2) ;                          % index of second points of gamma1 segements

hp = p(e2,:) - p(e1,:) ;
h = sqrt(sum(hp.^2,2)) ;                % length of the segments

[bu,bv,kbu,kbv] = buv(p,variables, refer) ;    % estimation of the boundary value at the nodes
bbu = h.*(bu(e1)+bu(e2))/2 ;                   % estimation of the boundary vector values for each segment
bbv = h.*(bv(e1)+bv(e2))/2 ;
kbbu = h.*(kbu(e1)+kbu(e2))/2 ;
kbbv = h.*(kbv(e1)+kbv(e2))/2 ;

% construction of the independent boundary vector
B = sparse(e1,1,bbu,2*np,1)+sparse(e2,1,bbu,2*np,1) + ...
    sparse(e1+np,1,bbv,2*np,1)+sparse(e2+np,1,bbv,2*np,1) ;

% construction of the linearly dependant (on Cu and Cv) boundary vector
Kb = sparse(e1,e1,kbbu,2*np,2*np)+sparse(e2,e2,kbbu,2*np,2*np) + ...
    sparse(e1+np,e1+np,kbbv,2*np,2*np)+sparse(e2+np,e2+np,kbbv,2*np,2*np) ;

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