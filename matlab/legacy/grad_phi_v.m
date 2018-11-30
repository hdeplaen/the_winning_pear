function [dphi1 , dphi2 , dphi3, T] = grad_phi_v(p,t,variables)
% Gradients of basis functions of Cv
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
%% RADIAL (supposedly not correct)
dphi1 = .5*[variables.Dvr*(-2*z32 + (r32.*z + r2.*z3 - r3.*z2)./r),  ...
    variables.Dvz*r32] ;
dphi2 = .5*[variables.Dvr*(2*z31 - (r31.*z + r1.*z3 - r3.*z1)./r), ...
    variables.Dvz*(-r31)] ;
dphi3 = .5*[variables.Dvr*(-2*z21 + (r21.*z + r1.*z2 - r2.*z1)./r),  ...
    variables.Dvz*(r21)] ;

%% EUCLIDEAN
% dphi1 = .5*[-variables.Dvr*z32, variables.Dvz*r32] ; 
% dphi2 = .5*[variables.Dvr*z31, -variables.Dvz*r31] ; 
% dphi3 = .5*[-z21*variables.Dvr, variables.Dvz*r21] ;

end