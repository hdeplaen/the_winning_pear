function [R] = stiff_matrix(p,t)
%STIFF_MATRIX Outputs stiff matrix
%   Henri De Plaen, KU Leuven
%
%  Input:
%         p : Node coordinates
%         t : Triangle vertices
%
%  Output:
%       R  : Stiffness matrix (symmetric positive semidefinite)
%            sparse np*np

np = size(p,1) ;

[dphi1, dphi2, dphi3] = grad_phi(p,t) ;

% cell-array of gradients
dphi = {dphi1 dphi2 dphi3} ;
R = sparse(np,np) ;

% under-diagonal entries
for i = 1:3
    for j=1:i-1
        R = R + sparse(t(:,i),t(:,j),sum(dphi{i}.*dphi{j},2),np,np) ;
    end
end

R = R + R.' ;

% diagonal entries
for i = 1:3
    R = R + sparse(t(:,i),t(:,i),sum(dphi{i}.*dphi{i},2),np,np) ;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dphi1 , dphi2 , dphi3] = grad_phi(p,t)
%GRADIENT PHI

r1 = p(t(:,1),1) ;
r2 = p(t(:,2),1) ;
r3 = p(t(:,3),1) ;

z1 = p(t(:,1),2) ;
z2 = p(t(:,2),2) ;
z3 = p(t(:,3),2) ;

r21 = r2-r1 ; z21 = z2-z1 ; 
r32 = r3-r2 ; z32 = z3-z2 ;
r31 = r3-r1 ; z31 = z3-z1 ;

% additions
z = (z1+z2+z3)/3 ;
r = (r1+r2+r3)/3 ;

% triangles area
T = (r21.*z31-z21.*r31)/2 ;

% gradients of basis functions
dphi1 = .5*[-2*z32 + (r32.*z + r2.*z3 - r3.*z2)./r,  r32] ;
dphi2 = .5*[ 2*z31 - (r31.*z + r1.*z3 - r3.*z1)./r, -r31] ;
dphi3 = .5*[-2*z21 + (r21.*z + r1.*z2 - r2.*z1)./r,  r21] ;

end

% r = full(stiff_matrix(p',t') );
