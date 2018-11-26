function [B,Kb] = boundary_vector_bis(p,e)
%BOUNDARY_VECTOR Outputs boundary vector
%   Henri De Plaen, KU Leuven
%
%  Input:
%         p : Node coordinates
%         e : Edge vector (bound1: index 1-5, bound2: index 6)
%         t : Triangle vertices
%
%  Output:
%       B  : Boundary vector
%            sparse np

%% FUNC BOUNDARY
np = size(p,1) ;

e6 = e(e(:,5)==6,:) ;
e1 = e6(:,1) ; 
e2 = e6(:,2) ;

hp = p(e2,:) - p(e1,:) ; 
h = sqrt(sum(hp.^2,2)) ;

[bu,bv,kbu,kbv] = buv(p) ;

bbu = h.*(bu(e1)+bu(e2))/2 ;
bbv = h.*(bv(e1)+bv(e2))/2 ;
kbbu = h.*(kbu(e1)+kbu(e2))/2 ;
kbbv = h.*(kbv(e1)+kbv(e2))/2 ;

B = sparse(e1,1,bbu,2*np,1)+sparse(e2,1,bbu,2*np,1) + ...
    sparse(e1+np,1,bbv,2*np,1)+sparse(e2+np,1,bbv,2*np,1) ;

Kb = sparse(e1,e1,kbbu,2*np,2*np)+sparse(e2,e2,kbbu,2*np,2*np) + ...
    sparse(e1+np,e1+np,kbbv,2*np,2*np)+sparse(e2+np,e2+np,kbbv,2*np,2*np) ;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bu,bv,kbu,kbv] = buv(p)

%% PARAMS (unused variables in comment)
% Dur = 2.8e-10 ;
% Duz = 1.1e-9 ;
% Dvr = 2.32e-9 ;
% Dvz = 6.97e-9 ;

Tcel = 25 ;
eta_u = 20.8/100 ;
eta_v = 0.04/100 ;

Rg = 8.314 ;
hu = 7e-7 ;
hv = 7.5e-7 ;
patm = 101300 ;
T = Tcel+273.15 ;
% Tref = 20 + 273.15 ;

Cuamb = patm*eta_u/Rg/T ;
Cvamb = patm*eta_v/Rg/T ;

% rq = 0.97 ;
% Kmfu = 0.1149 ;
% Kmv = 27.2438 ;
% Kmu = 0.4103 ;

% Eavmfvref = 56700 ;
% Vmfvref = 1.61e-4 ;
% Vmfv = Vmfvref*exp(Eavmfvref/Rg*(1/Tref-1/T)) ;

% Eavmuref = 80200 ;
% Vmuref = 2.39e-4 ;
% Vmu = Vmuref*exp(Eavmuref/Rg*(1/Tref-1/T)) ;

%% BOUNDARY VECTOR
r = p(:,1) ;
np = size(p,1) ;

% bu = r*-hu*Cuamb ;
% bv = r*hv*Cvamb ;
% 
% kbu = hu*r ;
% kbv = hv*r ;

bu = -hu*Cuamb.*ones(np,1) ;
bv = -hv*Cvamb.*ones(np,1) ;

kbu = hu.*ones(np,1) ;
kbv = hv.*ones(np,1) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%