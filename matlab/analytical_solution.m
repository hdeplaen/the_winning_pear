clear all ; clc ;

syms A B C R alpha beta gamma1 gamma2 r

f = A*cos(C*r) + B*sin(C*r) ;

df = diff(f,r) ;

d2f = 1/r^2*diff(r^2*df,r) ;

bnd = alpha*f+beta ;


r = 0 ; 
B = solve(subs(df)==0, B) ;

r = R ;
A = solve(subs(df)==subs(bnd), A) ;

syms r ;
%Rf = f ;
Rf = (beta*(2*sin(r) + r*cos(r)))/(r*sin(C*R) + alpha*r*cos(R)) ;
pretty(simplify(subs(d2f)))
C = solve(subs(d2f)==subs(Rf),C) ;

pretty(simplify(subs(f))) ;
