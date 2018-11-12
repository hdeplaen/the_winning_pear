% Author: Henri De Plaen, KU Leuven

clear all ; close all ; clc ;
syms a b c r1 r2 r3 z1 z2 z3 r z

for idx = 1:3

A = [r1 z1 1 ; r2 z2 1 ; r3 z3 1] ;
B = zeros(3,1) ; B(idx) = 1 ;
vsol = A\B ;

a = vsol(1) ; b = vsol(2) ; c = vsol(3) ;
T = (r1*z2 - r2*z1 - r1*z3 + r3*z1 + r2*z3 - r3*z2) ;
a = a*T ; b = b*T ; c = c*T ;

phi = a*r + b*z + c ;
sol = 1/r * diff(r*phi,r) ;
sol = simplify(expand(sol),10) ;
pretty(sol) ;

end