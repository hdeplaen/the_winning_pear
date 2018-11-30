syms Vmu Cu Kmu Cv Kmv rq Vmfv Kmfu ;

Ru = (C*beta*(2*sin(C*r) + C*r*cos(C*r)))/(C*r*sin(C*R) + alpha*r*cos(C*R)) ;
pretty(Ru) ;

Rv = (C*beta*(2*sin(C*r) + C*r*cos(C*r)))/(C*r*sin(C*R) + alpha*r*cos(C*R)) ;
pretty(Rv) ;

RudCu = simplify(diff(Ru,Cu))
RudCv = simplify(diff(Ru,Cv))

RvdCu = simplify(diff(Rv,Cu))
RvdCv = simplify(diff(Rv,Cv))

pretty(RudCu) ;
pretty(RudCv) ;
pretty(RvdCu) ;
pretty(RvdCv) ;