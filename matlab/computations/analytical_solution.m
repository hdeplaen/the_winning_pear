clear all ; clc ;

syms A B C Source R r hu Dur Cuamb r2

Cu = A*exp(-C*r)/r + B*exp(C*r)/C/r ;

dCu = diff(Cu,r) ;

d2Cu = diff(dCu,r) ;
d2Cur = 1/r*diff(r*dCu,r) ;

bnd = subs(hu/Dur*(Cu-Cuamb),r,R) ;

% r = 0 ;
% B = solve(subs(dCu)==0, B) ;
A = -B/C ;

% C = solve(subs(r^2*d2Cu+2*r*dCu) == subs(r^2*Source*Cu), C) ; % OK

r = R ; r2 = R ;
A = solve(subs(dCu)==subs(bnd), A) ;

syms r ;
subs(Cu) ;
