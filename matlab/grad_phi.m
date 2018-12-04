function varargout = grad_phi(varargin)
%Gradient of basis functions
%   [dphi1 , dphi2 , dphi3, T] = grad_phi(p,t,variables,refer,wvar)
%
%  Input:
%         p : Node coordinates
%         t : Triangle vertices
%         refer: 'cartesian' or 'radial'
%
%  Output:
%         dphi1 : gradient of first basis functions
%         dphi2 : gradient of second basis functions
%         dphi3 : gradient of thirs basis functions
%         T : area of the triangles
%
%   Henri De Plaen, KU Leuven

%% PRELIMINARIES
assert(nargin==3, 'Wrong number of input arguments') ;
p           = varargin{1} ;
t           = varargin{2} ;
refer       = varargin{3} ;

%% COMPUTE DISTANCES
r1 = p(t(:,1),1) ;
r2 = p(t(:,2),1) ;
r3 = p(t(:,3),1) ;

z1 = p(t(:,1),2) ;
z2 = p(t(:,2),2) ;
z3 = p(t(:,3),2) ;

r21 = r2-r1 ; z21 = z2-z1 ;
r32 = r3-r2 ; z32 = z3-z2 ;
r31 = r3-r1 ; z31 = z3-z1 ;

% additions
z = (z1+z2+z3)/3 ;
r = (r1+r2+r3)/3 ;

% triangles area
T = (r21.*z31-z21.*r31)/2 ;

%% RADIAL / EUCLID
switch refer
    case 'radial'
        %(supposedly not correct)
        dphi1 = .5*[(-2*z32 + (r32.*z + r2.*z3 - r3.*z2)./r), r32] ;
        dphi2 = .5*[(2*z31 - (r31.*z + r1.*z3 - r3.*z1)./r), -r31] ;
        dphi3 = .5*[(-2*z21 + (r21.*z + r1.*z2 - r2.*z1)./r),  r21] ;
        
    case 'cartesian'
        dphi1 = .5*[-z32, r32] ;
        dphi2 = .5*[z31, -r31] ;
        dphi3 = .5*[-z21, r21] ;
        
    otherwise ; error('Referential not recognized') ;
end

%% RETURN
assert(nargout==4, 'Wrong number of output arguments') ;

varargout{1} = dphi1 ;
varargout{2} = dphi2 ;
varargout{3} = dphi3 ;
varargout{4} = T ;

end