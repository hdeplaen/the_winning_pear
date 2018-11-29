function [B,stiff_bis] = boundary_vector_old(p,e,variables)
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
    [int1,int2,stiff11,stiff12,stiff21,stiff22] = int_sv_u(p1,p2,variables) ;
    
    B(l1,1) = B(l1,1) + int1 ;
    B(l2,1) = B(l2,1) + int2 ;
    
    stiff_bis(l1,l1) = stiff_bis(l1,l1) + stiff11 ;
    stiff_bis(l1,l2) = stiff_bis(l1,l2) + stiff12 ;
    stiff_bis(l2,l1) = stiff_bis(l2,l1) + stiff21 ;
    stiff_bis(l2,l2) = stiff_bis(l2,l2) + stiff22 ;
    
    %% Cv
    [int1,int2,stiff11,stiff12,stiff21,stiff22] = int_sv_v(p1,p2,variables) ;
    
    B(l1+np,1) = B(l1+np,1) + int1 ;
    B(l2+np,1) = B(l2+np,1) + int2 ;
    
    stiff_bis(l1+np,l1+np) = stiff_bis(l1+np,l1+np) + stiff11 ;
    stiff_bis(l1+np,l2+np) = stiff_bis(l1+np,l2+np) + stiff12 ;
    stiff_bis(l2+np,l1+np) = stiff_bis(l2+np,l1+np) + stiff21 ;
    stiff_bis(l2+np,l2+np) = stiff_bis(l2+np,l2+np) + stiff22 ;
end

end

% np = size(p,1) ; 
% ie1 = eneum(:,1) ; 
% ie2 = eneum(:,2) ;
% 
% % segments length
% xy = p(ie2,:)-p(ie1,:) ; 
% ee = sqrt(sum(xy.^2,2)) ;  
% gg = ee.*(g(ie1)+g(ie2))/2 ;    
% b = sparse(ie1,1,gg,np,1) + sparse(ie2,1,gg,np,1) ; 
% b = full(b) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [int1,int2,stiff11,stiff12,stiff21,stiff22] = int_sv_u(p1,p2,variables)
%% BOUNDARY VECTOR

h = sqrt(sum((p1-p2).^2)) ;
r1 = p1(1) ;
r2 = p2(1) ;

int1 = h/6*variables.hu*variables.Cuamb*(r1*2+r2) ;
int2 = h/6*variables.hu*variables.Cuamb*(r2*2+r1) ;

stiff11 = -h/6*variables.hu*2*r1 ;
stiff21 = -h/6*variables.hu*r2 ;
stiff12 = -h/6*variables.hu*2*r2 ;
stiff22 = -h/6*variables.hu*r1 ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [int1,int2,stiff11,stiff12,stiff21,stiff22] = int_sv_v(p1,p2,variables)
%% BOUNDARY VECTOR

h = sqrt(sum((p1-p2).^2)) ;
r1 = p1(1) ;
r2 = p2(1) ;

int1 = h/6*variables.hv*variables.Cvamb*(r1*2+r2) ;
int2 = h/6*variables.hv*variables.Cvamb*(r2*2+r1) ;

stiff11 = -h/6*variables.hv*2*r1 ;
stiff21 = -h/6*variables.hv*r2 ;
stiff12 = -h/6*variables.hv*2*r2 ;
stiff22 = -h/6*variables.hv*r1 ;
end