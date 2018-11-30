function [R] = stiff_matrix(p,t,variables,refer)
%STIFF_MATRIX Outputs stiff matrix
%   Henri De Plaen, KU Leuven
%
%  Input:
%         p : Node coordinates
%         t : Triangle vertices
%         variables : Problem parameters
%         refer : 'cartesian' or 'radial'
%
%  Output:
%       R  : Stiffness matrix (symmetric positive semidefinite)
%            sparse np*np

%% PRELIMINARIES
np = size(p,1) ; % number of nodes in total (boundary and non-boundary)
nt = size(t,1) ; % number of triangles

%% GRADIENTS OF BASIS FUNCTIONS
[dphi1, dphi2, dphi3, ~] = grad_phi(p,t,refer) ;   % gradients of basis functions

% cell-array of gradients
dphi = {dphi1 dphi2 dphi3} ;

Du = [variables.Dur*ones(nt,1) variables.Duz*ones(nt,1)];
Dv = [variables.Dvr*ones(nt,1) variables.Dvz*ones(nt,1)];

%% CONSTRUCTION OF THE STIFFNESS MATRIX
R = sparse(2*np,2*np) ; %prealloc

% under-diagonal entries
for i = 1:3
    for j=1:i-1
        R = R + sparse(t(:,i),t(:,j),sum(Du.*dphi{i}.*dphi{j},2),2*np,2*np) ;
        R = R + sparse(t(:,i)+np,t(:,j)+np,sum(Dv.*dphi{i}.*dphi{j},2),2*np,2*np) ;
    end
end

R = R + R.' ;

% diagonal entries
for i = 1:3
    R = R + sparse(t(:,i),t(:,i),sum(dphi{i}.*dphi{i},2),2*np,2*np) ;
    R = R + sparse(t(:,i)+np,t(:,i)+np,sum(dphi{i}.*dphi{i},2),2*np,2*np) ;
end

end