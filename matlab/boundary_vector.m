function [B,stiff_bis] = boundary_vector(p,e)
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

np = size(p,1) ;
ne = size(e,1) ;

% prealloc
B = sparse(2*np,1) ;
stiff_bis = sparse(2*np,2*np) ;

for idxe = find(e(5,:)~=6)
    loc_e = e(idxe,:) ;
    
    l1 = loc_e(1) ;
    l2 = loc_e(2) ;
    
    p1 = p(l1,:) ;
    p2 = p(l2,:) ;
    
    %% Cu
    [int1,int2,stiff11,stiff12,stiff21,stiff22] = int_sv_u(p1,p2) ;
    
    B(l1,1) = B(l1,1) + int1 ;
    B(l2,1) = B(l2,1) + int2 ;
    
    stiff_bis(l1,l1) = stiff_bis(l1,l1) + stiff11 ;
    stiff_bis(l1,l2) = stiff_bis(l1,l2) + stiff12 ;
    stiff_bis(l2,l1) = stiff_bis(l2,l1) + stiff21 ;
    stiff_bis(l2,l2) = stiff_bis(l2,l2) + stiff22 ;
    
    %% Cv
    [int1,int2,stiff11,stiff12,stiff21,stiff22] = int_sv_v(p1,p2) ;
    
    B(l1+np,1) = B(l1+np,1) + int1 ;
    B(l2+np,1) = B(l2+np,1) + int2 ;
    
    stiff_bis(l1+np,l1+np) = stiff_bis(l1+np,l1+np) + stiff11 ;
    stiff_bis(l1+np,l2+np) = stiff_bis(l1+np,l2+np) + stiff12 ;
    stiff_bis(l2+np,l1+np) = stiff_bis(l2+np,l1+np) + stiff21 ;
    stiff_bis(l2+np,l2+np) = stiff_bis(l2+np,l2+np) + stiff22 ;
end

end

np=size(p,1); ie1=eneum(:,1); ie2=eneum(:,2);
% segments length
xy=p(ie2,:)-p(ie1,:); ee=sqrt(sum(xy.^2,2));  
gg=ee.*(g(ie1)+g(ie2))/2;    
b=sparse(ie1,1,gg,np,1)+sparse(ie2,1,gg,np,1); 
b=full(b);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [int1,int2,stiff11,stiff12,stiff21,stiff22] = int_sv_u(p1,p2)

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

%% BOUNDARY VECTOR

h = sqrt(sum((p1-p2).^2)) ;
r1 = p1(1) ;
r2 = p2(1) ;

int1 = h/6*hu*Cuamb*(r1*2+r2) ;
int2 = h/6*hu*Cuamb*(r2*2+r1) ;

stiff11 = -h/6*hu*2*r1 ;
stiff21 = -h/6*hu*r2 ;
stiff12 = -h/6*hu*2*r2 ;
stiff22 = -h/6*hu*r1 ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [int1,int2,stiff11,stiff12,stiff21,stiff22] = int_sv_v(p1,p2)

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

%% BOUNDARY VECTOR

h = sqrt(sum((p1-p2).^2)) ;
r1 = p1(1) ;
r2 = p2(1) ;

int1 = h/6*hv*Cvamb*(r1*2+r2) ;
int2 = h/6*hv*Cvamb*(r2*2+r1) ;

stiff11 = -h/6*hv*2*r1 ;
stiff21 = -h/6*hv*r2 ;
stiff12 = -h/6*hv*2*r2 ;
stiff22 = -h/6*hv*r1 ;
end