function [dphi1 , dphi2 , dphi3, T] = grad_phi_u(p,t)
%% PARAMS (unused variables in comment)
Dur = 2.8e-10 ;
Duz = 1.1e-9 ;
Dvr = 2.32e-9 ;
Dvz = 6.97e-9 ;

Tcel = 25 ;
eta_u = 20.8/100 ;
eta_v = 0.04/100 ;

Rg = 8.314 ;
hu = 7e-7 ;
hv = 7.5e-7 ;
patm = 101300 ;
T = Tcel+273.15 ;
Tref = 20 + 273.15 ;

Cuamb = patm*eta_u/Rg/T ;
Cvamb = patm*eta_v/Rg/T ;

rq = 0.97 ;
Kmfu = 0.1149 ;
Kmv = 27.2438 ;
Kmu = 0.4103 ;

Eavmfvref = 56700 ;
Vmfvref = 1.61e-4 ;
Vmfv = Vmfvref*exp(Eavmfvref/Rg*(1/Tref-1/T)) ;

Eavmuref = 80200 ;
Vmuref = 2.39e-4 ;
Vmu = Vmuref*exp(Eavmuref/Rg*(1/Tref-1/T)) ;

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
dphi1 = .5*[Dur*(-2*z32 + (r32.*z + r2.*z3 - r3.*z2)./r),  Duz*r32] ;
dphi2 = .5*[Dur*(2*z31 - (r31.*z + r1.*z3 - r3.*z1)./r), Duz*(-r31)] ;
dphi3 = .5*[Dur*(-2*z21 + (r21.*z + r1.*z2 - r2.*z1)./r),  Duz*(r21)] ;

% dphi1 = .5*[-Dur*z32 Duz*r32] ; 
% dphi2 = .5*[Dur*z31 -Duz*r31] ; 
% dphi3 = .5*[-z21*Dur Duz*r21] ;

end

% r = full(stiff_matrix(p',t') );