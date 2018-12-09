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

B = zeros(2*np,1) ;                     % prealloc
Kb = zeros(2*np,2*np) ;                 % prealloc

g1 = e(:,5) ~= axis_e ;                 % index of gamma1 boundary segments (on the curve, not the symmetry axis)
l_g1 = length(g1) ;

e1 = e(g1,1) ;                          % index of first points of gamma1 segments
e2 = e(g1,2) ;                          % index of second points of gamma1 segements

hp = p(e2,:) - p(e1,:) ;
h = sqrt(sum(hp.^2,2)) ;                % length of the segments

for idx = 1:l_g1-1
    
    % BOUNDARY MATRIX
    Kb([g1(idx) g1(idx+1)],[g1(idx) g1(idx+1)]) = ...
        Kb([g1(idx) g1(idx+1)],[g1(idx) g1(idx+1)]) + ...
        kbu(r(idx:idx+1),l(idx), variables, refer) ;
    
    Kb([g1(idx) g1(idx+1)]+np,[g1(idx) g1(idx+1)]+np) = ...
        Kb([g1(idx) g1(idx+1)]+np,[g1(idx) g1(idx+1)]+np) + ...
        kbv(r(idx:idx+1),l(idx), variables, refer) ;
    
    % BOUNDARY VECTOR
    B([g1(idx) g1(idx+1)]) = B([g1(idx) g1(idx+1)]) + ...
        bu(e1(idx:idx+1), h(idx), variables, refer) ;
    
    B([g1(idx) g1(idx+1)]+np) = B([g1(idx) g1(idx+1)]+np) + ...
        bv(e1(idx:idx+1), h(idx), variables, refer) ;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bu_out = bu(r,h,variables, refer)

switch refer
    case 'radial'
        bu_out = (h/6)*[ 2*r(1) + r(2) ; r(1) + 2*r(2) ] ;
    case 'cartesian'
        error('Not computed for the moment') ;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%