function [R] = stiff_matrix(p,t,variables,refer)
%STIFF_MATRIX Outputs stiff matrix
%   Henri De Plaen, KU Leuven
%
%  Input:
%         p : Node coordinates
%         t : Triangle vertices
%         variables : Problem parameters
%         refer : 'euclid' or 'radial'
%
%  Output:
%       R  : Stiffness matrix (symmetric positive semidefinite)
%            sparse np*np

%% PRELIMINARIES
np = size(p,1) ;                                            % number of nodes in total (boundary and non-boundary)

%% GRADIENTS OF BASIS FUNCTIONS
[dphi1u, dphi2u, dphi3u, ~] = grad_phi(p,t,variables,refer,'u') ;   % gradients of basis functions of Cu
[dphi1v, dphi2v, dphi3v, ~] = grad_phi(p,t,variables,refer,'v') ;   % gradients of basis functions of Cv

% cell-array of gradients
dphiu = {dphi1u dphi2u dphi3u} ;
dphiv = {dphi1v dphi2v dphi3v} ;

%% CONSTRUCTION OF THE STIFFNESS MATRIX
R = sparse(2*np,2*np) ; %prealloc

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