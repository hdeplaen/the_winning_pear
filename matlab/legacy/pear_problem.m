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
max_it = 1e+3 ;
tol_min = 1e-10 ;

%% TYPE OF PROBLEM
variables = generate_variables(type) ;                      % problem variables
[p,e,t,model,axis_s] = generate_mesh(hmax,'test',0,0) ;     % mesh parameters

K = stiff_matrix_new(p,t,variables,refer) ;                 % compute stiff matrix
[B,Rb] = boundary_vector_new(p,e,variables,axis_s,refer) ;  % compute boundary

np = size(p,1) ;                                            % number of points
sol_buff = abs(.1*randn(2*np,1)) ;                          % generate initial guess

%% SOLVE
%FIRST ITERATION
[F,dF] = fun_vector_new(p,t,sol_buff,variables,refer) ;     % compute source terms
dF = 1/regv_init*dF ;                                       % ajust step size
sol_buff = (K+Rb-dF)\(F-dF*sol_buff+B) ;                    % compute solution

%NEWTON-RAPHSON EXPLICIT SOVER
h = waitbar(0,'Initializing') ; % creating a waitbar

for it = 1:max_it
    [F,dF] = fun_vector_new(p,t,sol_buff,variables,refer) ;     % compute source terms
    dF = 1/regv*dF ;                                            % ajust step size
    sol = (K+Rb-dF)\(F-dF*sol_buff+B) ;                         % compute solution
    
    % verify tolerance
    tol = norm(sol-sol_buff) ;
    if tol < tol_min ; break ; end
    
    % update buffer
    sol_buff = sol ;
    
    % update waitbar
    if mod(it,300)==0
        waitbar(it/max_it,h,['tol = ' num2str(tol) newline ...
            'iteration = ' num2str(it) ' / ' num2str(max_it)]) ;
    end
end

delete(h) ;

sol_u = sol(1:np) ;
sol_v = sol(np+1:end) ;

%% VIEW SOL
if plot_sol    
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
varargout{1} = sol_u ;
varargout{2} = sol_v ;

if strcmp(type,'test')
    warning('not implemented now') ;
end

end