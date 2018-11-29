function [Ru,Rv,RudCu,RudCv,RvdCu,RvdCv] = fun_diff_test(sol, variables)
%FUN_DIFF Henri De Plaen

n_sol = size(sol,1) ;

assert(mod(n_sol,2)==0) ;
Cu = sol(1:n_sol/2,:) ;
Cv = sol(n_sol/2+1:end,:) ;

%% FUN
Ru = -variables.Vmu.*Cu./((variables.Kmu+Cu).*(1+Cv./variables.Kmv)) ;
Rv = -variables.rq.*Ru + variables.Vmfv./(1+Cu./variables.Kmfu) ;

RudCu = -(variables.Kmu.*variables.Kmv.*variables.Vmu)./ ...
    ((Cu + variables.Kmu).^2.*(Cv + variables.Kmv)) ;
RudCv = (Cu.*variables.Kmv.*variables.Vmu)./ ...
    ((Cu + variables.Kmu).*(Cv + variables.Kmv).^2) ; 
RvdCu = (variables.Kmv.*variables.Vmu.*variables.rq)./ ...
    ((Cu + variables.Kmu).*(Cv + variables.Kmv)) - ...
    (variables.Kmfu.*variables.Vmfv)./(Cu + variables.Kmfu).^2 - ...
    (Cu.*variables.Kmv.*variables.Vmu.*variables.rq)./ ...
    ((Cu + variables.Kmu).^2.*(Cv + variables.Kmv)) ;
RvdCv = -(Cu.*variables.Kmv.*variables.Vmu.*variables.rq)./ ...
    ((Cu + variables.Kmu).*(Cv + variables.Kmv).^2) ;

end

