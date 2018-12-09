function varargout = generate_mesh_native(varargin)
%Generates mesh for the pear problem
%   [p,e,t] = generate_mesh(hmax,type,export_csv,plot_mesh)
%
% INPUT
%   hmax: maximum edge length (controls the mesh concentration)
%   type: 'pear' or 'test'
%   export_csv: create csv (true or false)
%   plot_mesh: plot mesh (true or false)
% OUTPUT
%   p: mesh points (2*np matrix)
%   e: mesh edges (3*ne matrix)
%   t: mesh triangles (2*nt matrix)
%
%Author: Henri De Plaen
%KU Leuven
%Project WIT: the winning pear
%Date: Nov 2018

%% PRELIMINARIES
assert(nargin==4, 'Wrong number of input arguments (4)' ) ;
hmax = varargin{1} ;
type = varargin{2} ;
export_csv = varargin{3} ;
plot_mesh = varargin{4} ;

%% MESH
model = createpde(1) ;

% create geometry of good model
switch type
    case 'test'
        C1 = [1,0,1,1]';
        R1 = [3,4,-1,-1,0,0,2,0,0,2]' ;
        C1 = [C1;zeros(length(R1) - length(C1),1)];
        gm = [R1,C1];
        sf = 'C1-R1';
        ns = char('R1','C1');
        ns = ns';
        g = decsg(gm,sf,ns);
        geometryFromEdges(model,g) ;
        
        axis_s = 1 ;
    case 'pear'
        geometryFromEdges(model,@pear_geometry) ;
        axis_s = 6 ;
    otherwise
        error('type not recongized') ;
end



% generate mesh
mesh = generateMesh(model,'GeometricOrder','linear', 'Hgrad', 1.9, 'Hmax', hmax) ;
[p,e_old,t] = meshToPet(mesh);

e = e_old(:,1:2) ;
e(:,3) = 1 ;
idx_s = e_old(:,5) == axis_s ;
e(idx_s,3) = 2 ; 

p = p/50 ;

%% PLOT
if plot_mesh
    pdemesh(model,'NodeLabels','off') ; %,'ElementLabels','on') ;
    title('Pear mesh') ;
    xlabel('x') ; ylabel('y') ;
    axis([-.1 0.35 -.1 1.1]) ;
end

%% EXPORT CSV
if export_csv
    % format needed for C++ program (to be completed later)
    node_csv = t(1:3,:)'-1 ;
    xp_csv = p(1,:)' ;
    yp_csv = p(2,:)' ;
    
    csvwrite('node.csv',node_csv) ;
    csvwrite('xp.csv',xp_csv) ;
    csvwrite('yp.csv',yp_csv) ;
end

%% RETURN
assert(nargout==3, 'Wrong number of input parameters (3)') ;

varargout{1} = p' ;
varargout{2} = e' ;
varargout{3} = t' ;

end