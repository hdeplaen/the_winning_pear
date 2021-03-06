function [K] = stiff_matrix(p,t,variables,~)
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

%% CONSTKuCTION OF THE STIFFNESS MATRIX
Ku = zeros(np,np) ; % prealloc
Kv = zeros(np,np) ; % prealloc

for m = 1:length(t)
    % nodes of the mesh element (triangle)
    t_loc = t(m,1:3) ;
    r = p(t_loc,1) ;
    z = p(t_loc,2) ;
    
    Ku(t_loc,t_loc) = Ku(t_loc,t_loc) + ...
        stiff_block(r,z,variables.Dur,variables.Duz);
    Kv(t_loc,t_loc) = Kv(t_loc,t_loc) + ...
        stiff_block(r,z,variables.Dvr,variables.Dvz);
end

K = [sparse(Ku) sparse(np,np) ; sparse(np,np) sparse(Kv)] ;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ K_block ] = stiff_block(xp,yp, Dr, Dz)
[T,~,dphi2,dphi3] = grad_phi(xp,yp) ;
K_block = zeros(3,3) ;

for idx1 = 1:3
   for idx2 = 1:3       
       K_block(idx1,idx2) = (sum(xp)).* ...
           (Dr*dphi2(idx1)*dphi2(idx2) + ...
           Dz*dphi3(idx1)*dphi3(idx2))./(12*T) ;
   end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T,dphi1,dphi2,dphi3] = grad_phi(xp,yp)
    dphi1 = [] ;
    dphi2 = circshift(yp,2)-circshift(yp,1) ;
    dphi3 = circshift(xp,1)-circshift(xp,2) ;
    
    % mesh element area
    T = xp(2)*yp(3) + xp(1)*yp(2) + ...
        xp(3)*yp(1) - xp(2)*yp(1) - ...
        xp(1)*yp(3) - xp(3)*yp(2) ; 
end