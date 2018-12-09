syms Vmu Cu Kmu Cv Kmv rq Vmfv Kmfu ;

Ru = Vmu*Cu/((Kmu+Cu)*(1+Cv/Kmv)) ;
pretty(Ru) ;

Rv = rq*Ru + Vmfv/(1+Cu/Kmfu) ;
pretty(Rv) ;

RudCu = simplify(diff(Ru,Cu))
RudCv = simplify(diff(Ru,Cv))

RvdCu = simplify(diff(Rv,Cu))
RvdCv = simplify(diff(Rv,Cv))

pretty(RudCu) ;
pretty(RudCv) ;
pretty(RvdCu) ;
pretty(RvdCv) ;