function [dphi1 , dphi2 , dphi3, T] = grad_phi_u(p,t,variables)
assert(nargin==3, 'Wrong number of input arguments') ;

%% GRADIENT PHI

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
dphi1 = .5*[variables.Dur*(-2*z32 + (r32.*z + r2.*z3 - r3.*z2)./r), ...
    variables.Duz*r32] ;
dphi2 = .5*[variables.Dur*(2*z31 - (r31.*z + r1.*z3 - r3.*z1)./r), ...
    variables.Duz*(-r31)] ;
dphi3 = .5*[variables.Dur*(-2*z21 + (r21.*z + r1.*z2 - r2.*z1)./r),  ...
    variables.Duz*(r21)] ;

% dphi1 = .5*[-variables.Dur*z32, variables.Duz*r32] ; 
% dphi2 = .5*[variables.Dur*z31, -variables.Duz*r31] ; 
% dphi3 = .5*[-z21*variables.Dur, variables.Duz*r21] ;

assert(nargout==4, 'Wrong number of output arguments') ;

end

% r = full(stiff_matrix(p',t') );