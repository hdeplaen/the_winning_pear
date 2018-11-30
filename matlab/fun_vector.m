function [F,dF] = fun_vector(p,t,sol,variables,refer)
%Outputs function vector for linearization F(x) + dF(x)(x1-x0)
%   [F,dF] = fun_vector(p,t,sol,variables)
%
%  Input:
%         p : Node coordinates
%         sol : points at which F(sol) + dF(sol) should be evaluated
%         variables : problem parameters
%         refer : 'euclid' or 'radial'
%
%  Output:
%       F  : Function vector
%            sparse (2*np,1)
%       dF : Jacobian of the function vector
%            sparse (2*np,2*np)
%
% Henri De Plaen, KU Leuven

np = size(p,1) ;

% prealloc
F = sparse(2*np,1) ;
dF = sparse(2*np,2*np) ;

% comuting the source term and jacobian at each node
[Ru,Rv,RudCu,RudCv,RvdCu,RvdCv] = fun_diff(sol,variables) ;

% evaluation of the source term on the triangles
fu = (Ru(t(:,1))+Ru(t(:,2))+Ru(t(:,3)))/3 ;
fv = (Rv(t(:,1))+Rv(t(:,2))+Rv(t(:,3)))/3 ;

% evaluation of the jacobian of the source term on the triangles
fudu = (RudCu(t(:,1))+RudCu(t(:,2))+RudCu(t(:,3)))/3 ;
fudv = (RudCv(t(:,1))+RudCv(t(:,2))+RudCv(t(:,3)))/3 ;
fvdu = (RvdCu(t(:,1))+RvdCu(t(:,2))+RvdCu(t(:,3)))/3 ;
fvdv = (RvdCv(t(:,1))+RvdCv(t(:,2))+RvdCv(t(:,3)))/3 ;

% area of the triangles
[~ , ~ , ~, T] = grad_phi(p,t,variables,refer,'u') ;

% division by the triangles area
fu = fu.*T/3 ;
fv = fv.*T/3 ;
fudu = fudu.*T/3 ;
fudv = fudv.*T/3 ;
fvdu = fvdu.*T/3 ;
fvdv = fvdv.*T/3 ;

% construction of the source term F and dF by adding each triangle
% component at each corresponding node
for idx=1:3
    F = F + sparse(t(:,idx),1,fu,2*np,1) ;
    F = F + sparse(t(:,idx)+np,1,fv,2*np,1) ;
    
    dF = dF + sparse(t(:,idx),t(:,idx),fudu,2*np,2*np) ;
    dF = dF + sparse(t(:,idx),t(:,idx)+np,fudv,2*np,2*np) ;
    dF = dF + sparse(t(:,idx)+np,t(:,idx),fvdu,2*np,2*np) ;
    dF = dF + sparse(t(:,idx)+np,t(:,idx)+np,fvdv,2*np,2*np) ;
end

end