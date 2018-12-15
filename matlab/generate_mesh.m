function varargout = generate_mesh(varargin)
%Generates mesh for the pear problem
%   [p,e,t] = generate_mesh(hmax,type,export_csv,plot_mesh)
%
% INPUT
%   hmax        : maximum edge length (controls the mesh concentration)
%   type        : 'pear' or 'test'
%   export_csv  : create csv (true or false)
%   plot_mesh   : plot mesh (true or false)
% OUTPUT
%   p           : mesh points (2*np matrix)
%   e           : mesh edges (3*ne matrix) with e(3,:) index of boundary (1 for
%                 exterior, 2 for axis)
%   t           : mesh triangles (4*nt matrix)
%
%Author: Henri De Plaen
%KU Leuven
%Project WIT: the winning pear
%Date: Dec 2018

%% PRELIMINARIES
assert(nargin==4, 'Wrong number of input arguments (4)' ) ;
hmax = varargin{1} ;
type = varargin{2} ;
export_csv = varargin{3} ;
plot_mesh = varargin{4} ;

%% BOUNDARY
switch type
    case 'test'
        ep = half_circle(round(1/hmax)) ;
        
    case 'pear'
        ep = eval_geometry(@pear_geometry,round(1/hmax)) ;
    otherwise
        error('type not recongized') ;
end

%% EVULATE MESH
figure ;
[p,t]=distmesh_2d(@dpoly,@huniform,hmax,[-1,-1; 1,2],1e+3,ep,ep);
e_temp = boundary(p(:,1),p(:,2)) ;

e(1,:) = e_temp(1:end-1)' ;
e(2,:) = e_temp(2:end)' ;
e(3,:) = 1 ;

log_temp = p(e_temp,1) < 1e-3 ;
idx_g2 = or(log_temp(1:end-1),log_temp(2:end)) ;
e(3,idx_g2) = 2 ;


%% PLOT
if plot_mesh
    pdemesh(model,'NodeLabels','on') ; %,'ElementLabels','on') ;
    title('Pear mesh') ;
    xlabel('x') ; ylabel('y') ;
    axis([-.1 0.35 -.1 1.1]) ;
end

%% EXPORT CSV
if export_csv
    % format needed for C++ program (to be completed later)
    node_csv = t(:,1:3)'-1 ;
    points_csv = p(:,1:2)' ;
    boundary_csv = e ;
    
    csvwrite('exports/node.csv',node_csv) ;
    csvwrite('exports/points.csv',points_csv) ;
    csvwrite('exports/boundary.csv',boundary_csv) ;
end

%% RETURN
assert(nargout==3, 'Wrong number of output arguments (3)') ;

varargout{1} = p ;
varargout{2} = e' ;
varargout{3} = t ;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ep = eval_geometry(geom_function, prec)

ne = feval(geom_function) ;                     % number of edges
prec = round(prec/ne/12) ;

ep = zeros(ne*prec,2) ;                      % prealloc
s = linspace(0,1,prec) ;                       % arc range

for idx = 1:ne
    [x,y] = feval(geom_function,idx,s) ;
    ep( (idx-1)*prec+1 : (idx)*prec, 1) = x(:) ;
    ep( (idx-1)*prec+1 : (idx)*prec, 2) = y(:) ;
end

ep = ep/20 ;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ep = half_circle(prec)

n_points_circle  = 32 ;
n_points_quarter = n_points_circle/2 ;
n_points_total   = 50 ;

r = zeros(n_points_total+1,1) ;         % prealloc
z = zeros(n_points_total+1,1) ;         % prealloc

for idx = 2:n_points_total+1
    if idx <= n_points_circle + 1
        r(idx) = sin((idx-1) * pi/2/n_points_quarter) ;
        z(idx) = cos((idx-1) * pi/2/n_points_quarter) ;
     elseif idx > n_points_circle + 1
        r(idx) = 0;
        z(idx) = -1 + 2 * (idx-n_points_circle-1)/ ...
            (n_points_total-n_points_circle) ;
    end
end

ep(:,1) = r(:) ;
ep(:,2) = z(:) ;
ep(1,:) = [0 1] ;

ep = ep/20 ;             % resizing

end