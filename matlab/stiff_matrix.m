function [R] = stiff_matrix(p,t,variables,~)
%STIFF_MATRIX Outputs stiff matrix
%   Henri De Plaen, KU Leuven
%
%  Input:
%         p : Node coordinates
%         t : Triangle vertices
%         variables : Problem parameters
%         refer : 'cartesian' or 'radial' (legacy)
%
%  Output:
%       R  : Stiffness matrix (symmetric positive semidefinite)
%            sparse np*np

%% PRELIMINARIES
np = size(p,1) ; % number of nodes in total (boundary and non-boundary)
nt = size(t,1) ; % number of triangles

%% GRADIENTS OF BASIS FUNCTIONS
% Du = [variables.Dur*ones(nt,1) variables.Duz*ones(nt,1)];
% Dv = [variables.Dvr*ones(nt,1) variables.Dvz*ones(nt,1)];

%% CONSTRUCTION OF THE STIFFNESS MATRIX
Ru = zeros(np,np) ; % prealloc
Rv = zeros(np,np) ; % prealloc

for m = 1:length(t)
    % nodes of the mesh element (triangle)
    t_loc = t(m,1:3) ;
    r = p(t_loc,1) ;
    z = p(t_loc,2) ;
    
    Ru(t_loc,t_loc) = Ru(t_loc,t_loc) + ...
        stiff_block(r,z,variables.Dur,variables.Duz);
    Rv(t_loc,t_loc) = Rv(t_loc,t_loc) + ...
        stiff_block(r,z,variables.Dvr,variables.Dvz);
end

R = [sparse(Ru) sparse(np,np) ; sparse(np,np) sparse(Rv)] ;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ R_block ] = stiff_block(r,z, Dr, Dz)
[T,~,dphi1,dphi2] = grad_phi(r,z) ;
R_block = zeros(3,3) ;

for i = 1:3
   for j = 1:3       
       R_block(i,j) = (r(1) + r(2) + r(3))./(12*T).* ...
           (Dr*dphi1(i)*dphi1(j) + Dz*dphi2(i)*dphi2(j)) ;
   end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T,dphi1,dphi2,dphi3] = grad_phi(r,z)
    dphi1 = [] ;
    dphi2 = circshift(z,2)-circshift(z,1) ;
    dphi3 = circshift(r,1)-circshift(r,2) ;
    
    % mesh element area
    T = r(2)*z(3)+r(1)*z(2) + ...
        r(3)*z(1) - r(2)*z(1) - ...
        r(1)*z(3) - r(3)*z(2) ; 
end