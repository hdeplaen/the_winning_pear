function [R] = stiff_matrix(p,t,e)
%STIFF_MATRIX Outputs stiff matrix
%   Henri De Plaen, KU Leuven
%
%  Input:
%         p : Node coordinates
%         t : Triangle vertices
%
%  Output:
%       R  : Stiffness matrix (symmetric positive semidefinite)
%            sparse np*np

np = size(p,1) ;

[dphi1u, dphi2u, dphi3u] = grad_phi_u(p,t) ;
[dphi1v, dphi2v, dphi3v] = grad_phi_v(p,t) ;

% cell-array of gradients
dphiu = {dphi1u dphi2u dphi3u} ;
dphiv = {dphi1v dphi2v dphi3v} ;
R = sparse(2*np,2*np) ;

% under-diagonal entries
for i = 1:3
    for j=1:i-1
        R = R + sparse(t(:,i),t(:,j),sum(dphiu{i}.*dphiu{j},2),2*np,2*np) ;
        R = R + sparse(t(:,i)+np,t(:,j)+np,sum(dphiv{i}.*dphiv{j},2),2*np,2*np) ;
    end
end

R = R + R.' ;

% diagonal entries
for i = 1:3
    R = R + sparse(t(:,i),t(:,i),sum(dphiu{i}.*dphiu{i},2),2*np,2*np) ;
    R = R + sparse(t(:,i)+np,t(:,i)+np,sum(dphiv{i}.*dphiv{i},2),2*np,2*np) ;
end

end
