function varargout = pear_problem(varargin)
%Generates mesh for the pear problem
%   [p,e,t] = pear_problem(hmax,type,plot_sol)
%
% INPUT
%   hmax: maximum edge length (controls the mesh concentration)
%   type: 'pear' or 'test'
%   plot_sol: plot solution (true or false)
%   refer: 'cartesian' or 'radial'
% OUTPUT
%   sol: solution [Cu; Cv]
%   err: L2 error in case of the test
%
%Author: Henri De Plaen
%KU Leuven
%Project WIT: the winning pear
%Date: Nov 2018

%% PRELIMINARIES
assert(nargin==4, 'Wrong number of input arguments (4)' ) ;
hmax        = varargin{1} ;
type        = varargin{2} ;
plot_sol    = varargin{3} ;
refer       = varargin{4} ;

% set seed for reproducability of experiments
% rng(0681349)

% step sizes
regv_init = 1e+1 ;
regv = 1e+0 ;

% non-linear slover parameters
max_it = 1e+4 ;
tol_min = 1e-10 ;

%% TYPE OF PROBLEM
variables = generate_variables(type) ;                  % problem variables
[p,e,t,model,axis_s] = generate_mesh(hmax,type,0,0) ;   % mesh parameters

K = stiff_matrix(p,t,variables,refer) ;                 % compute stiff matrix
e_curve = e(e(:,5)~=axis_s,:) ;                         % boundary not on the axis (Gamma1)
ibcd = union(e_curve(:,1),e_curve(:,2)) ;               % index of boundary elements

[B,Rb] = boundary_vector(p,e,variables,axis_s,refer) ;  % compute boundary

np = size(p,1) ;                                        % number of points
sol_buff = abs(.1*randn(2*np,1)) ;                     % generate initial guess

%% SOLVE
%FIRST ITERATION
[F,dF] = fun_vector(p,t,sol_buff,variables,refer) ;     % compute source terms
F([ibcd ibcd+np]) = 0 ;                                 % no source term on the boundary (F)
dF([ibcd ibcd+np],[ibcd ibcd+np]) = 0 ;                 % no source term on the boundary (dF)
dF = 1/regv_init*dF ;                                   % ajust step size
sol = (K-Rb-dF)\(F-dF*sol_buff+B) ;                     % compute solution

%NEWTON-RAPHSON EXPLICIT SOVER
h = waitbar(0,'Initializing') ; % creating a waitbar

for it = 1:max_it
    [F,dF] = fun_vector(p,t,sol_buff,variables,refer) ;    % compute source terms
    F([ibcd ibcd+np]) = 0 ;                                % no source term on the boundary (F)
    dF([ibcd ibcd+np],[ibcd ibcd+np]) = 0 ;                % no source term on the boundary (dF)
    dF = 1/regv*dF ;                                       % ajust step size
    sol = (K-Rb-dF)\(F-dF*sol_buff+B) ;                    % compute solution
    
    % verify tolerance
    tol = sqrt(mean((sol-sol_buff).^2)) ;
    if abs(tol) < tol_min ; break ; end
    
    % update buffer
    sol_buff = sol ;
    
    % update waitbar
    if mod(it,300)==0
        waitbar(it/max_it,h,['tol = ' num2str(tol) newline ...
            'iteration = ' num2str(it) ' / ' num2str(max_it)]) ;
    end
end

delete(h) ;

%% VIEW SOL
if plot_sol
    sol_u = sol(1:np) ;
    sol_v = sol(np+1:end) ;
    
    figure ; hold on ;
    
    subplot(1,2,1) ;
    title('Oxygen concentration') ;
    pdeplot(model,'XYData',sol_u(:)) ;
    xlabel('r') ; ylabel('z') ;
    colormap('parula') ;
    
    subplot(1,2,2) ;
    title('Carbon dioxide concentration') ;
    pdeplot(model,'XYData',sol_v(:)) ;
    xlabel('r') ; ylabel('z') ;
    colormap('parula') ;
end

%% RETURN (not ready for the moment)
%assert(nargout>2, 'Too much output arguments (>2)') ;
l2_err = 0 ; % not defined for the moment
varargout{1} = sol ;
varargout{2} = l2_err ;
end