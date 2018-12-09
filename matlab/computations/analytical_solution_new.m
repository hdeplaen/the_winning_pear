clear all ; clc ;

syms r A B hu Dur Cuamb Source R

Cu = A*exp(-Source*r)/r + B*exp(Source*r)/Source/r ;

dCu = diff(Cu,r) ;

bnd = hu/Dur*(Cu-Cuamb) ;

r = R ;
A = solve(subs(dCu)==subs(bnd), A) ;

r = 0 ; 
B = solve(subs(dCu)==0, B) ;

syms r ;
%Rf = f ;
Rf = (beta*(2*sin(r) + r*cos(r)))/(r*sin(C*R) + alpha*r*cos(R)) ;
pretty(simplify(subs(d2f)))
% C = solve(subs(d2f)==subs(Rf),C) ;
% 
% pretty(simplify(subs(f))) ;
