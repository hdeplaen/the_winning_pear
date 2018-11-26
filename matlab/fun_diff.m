function [Ru,Rv,RudCu,RudCv,RvdCu,RvdCv] = fun_diff(sol)
%FUN_DIFF Henri De Plaen

n_sol = size(sol,1) ;

assert(mod(n_sol,2)==0) ;
Cu = sol(1:n_sol/2,:) ;
Cv = sol(n_sol/2+1:end,:) ;

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

%% FUN
Ru = -Vmu.*Cu./((Kmu+Cu).*(1+Cv./Kmv)) ;
Rv = -rq.*Ru + Vmfv./(1+Cu./Kmfu) ;

% RudCu = (Kmu.*Kmv.*Vmu)./((Cu + Kmu)..^2.*(Cv + Kmv)) ;  
% RudCv = -(Cu.*Kmv.*Vmu)./((Cu + Kmu).*(Cv + Kmv)..^2) ;
% RvdCu = (Kmv.*Vmu.*rq)./((Cu + Kmu).*(Cv + Kmv)) - ...
%     (Kmfu.*Vmfv)./(Cu + Kmfu)..^2 - ...
%     (Cu.*Kmv.*Vmu.*rq)./((Cu + Kmu)..^2.*(Cv + Kmv)) ; 
% RvdCv = -(Cu.*Kmv.*Vmu.*rq)./((Cu + Kmu).*(Cv + Kmv)..^2) ;

RudCu = -(Kmu.*Kmv.*Vmu)./((Cu + Kmu).^2.*(Cv + Kmv)) ;
RudCv = (Cu.*Kmv.*Vmu)./((Cu + Kmu).*(Cv + Kmv).^2) ; 
RvdCu = (Kmv.*Vmu.*rq)./((Cu + Kmu).*(Cv + Kmv)) - ...
    (Kmfu.*Vmfv)./(Cu + Kmfu).^2 - ...
    (Cu.*Kmv.*Vmu.*rq)./((Cu + Kmu).^2.*(Cv + Kmv)) ;
RvdCv = -(Cu.*Kmv.*Vmu.*rq)./((Cu + Kmu).*(Cv + Kmv).^2) ;

end

