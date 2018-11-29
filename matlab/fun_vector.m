function [F,dF] = fun_vector(p,t,sol,variables)
%FUN VECTOR Outputs function vector
%   Henri De Plaen, KU Leuven
%
%  Input:
%         p : Node coordinates
%
%  Output:
%       F  : Function vector
%            sparse np

np = size(p,1) ;

% prealloc
F = sparse(2*np,1) ;
dF = sparse(2*np,2*np) ;

[Ru,Rv,RudCu,RudCv,RvdCu,RvdCv] = fun_diff(sol,variables) ;
[~ , ~ , ~, T] = grad_phi_u(p,t,variables) ;

fu = (Ru(t(:,1))+Ru(t(:,2))+Ru(t(:,3)))/3 ;
fv = (Rv(t(:,1))+Rv(t(:,2))+Rv(t(:,3)))/3 ;
fudu = (RudCu(t(:,1))+RudCu(t(:,2))+RudCu(t(:,3)))/3 ;
fudv = (RudCv(t(:,1))+RudCv(t(:,2))+RudCv(t(:,3)))/3 ;
fvdu = (RvdCu(t(:,1))+RvdCu(t(:,2))+RvdCu(t(:,3)))/3 ;
fvdv = (RvdCv(t(:,1))+RvdCv(t(:,2))+RvdCv(t(:,3)))/3 ;

fu = fu.*T/3 ;
fv = fv.*T/3 ;
fudu = fudu.*T/3 ;
fudv = fudv.*T/3 ;
fvdu = fvdu.*T/3 ;
fvdv = fvdv.*T/3 ;

for idx=1:3
    F = F + sparse(t(:,idx),1,fu,2*np,1) ;
    F = F + sparse(t(:,idx)+np,1,fv,2*np,1) ;
    
    dF = dF + sparse(t(:,idx),t(:,idx),fudu,2*np,2*np) ;
    dF = dF + sparse(t(:,idx),t(:,idx)+np,fudv,2*np,2*np) ;
    dF = dF + sparse(t(:,idx)+np,t(:,idx),fvdu,2*np,2*np) ;
    dF = dF + sparse(t(:,idx)+np,t(:,idx)+np,fvdv,2*np,2*np) ;
end

end