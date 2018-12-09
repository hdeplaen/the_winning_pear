function varargout = pear_problem(varargin)
%Generates mesh for the pear problem
%   [Cu, Cv] = pear_problem(hmax, type, plot_sol, refer)
%
% INPUTS
%   hmax        : maximum edge length (controls the mesh concentration)
%   type        : 'pear' or 'test'
%   plot_sol    : plot solution (true or false)
%   refer       : 'cartesian' or 'radial'
%   Tcel        : parameter 1
%   eta_u       : parameter 2
%   eta_v       : parameter 3
%   casus       : name of the studied case
%
% OUTPUTS
%   Cu          : solution for Cu
%   Cv          : solution for Cv
%
%Author: Henri De Plaen
%KU Leuven
%Project WIT: the winning pear
%Date: Dec 2018

%% PRELIMINARIES
assert(nargin==8, 'Wrong number of input arguments (8)' ) ;
hmax        = varargin{1} ;
type        = varargin{2} ;
plot_sol    = varargin{3} ;
refer       = varargin{4} ;
Tcel        = varargin{5} ;
eta_u       = varargin{6} ;
eta_v       = varargin{7} ;
casus       = varargin{8} ;

% non-linear slover parameters
max_it = 1e+3 ;
tol_min = 1e-10 ;

%% TYPE OF PROBLEM
variables = generate_variables('pear',Tcel, eta_u, eta_v) ;                      % problem variables
[p,e,t] = generate_mesh(hmax,type,0,0) ;                      % mesh parameters

K = stiff_matrix(p,t,variables,refer) ;                       % compute stiff matrix
[Bu,Bv,Kbu,Kbv] = boundary_vector(p,e,variables,refer) ;      % compute boundary

np = size(p,1) ;                                              % number of points

%% SOLVE
%LINEAR SOLUTION
Ku = K(1:np,1:np) ;
Kv = K(np+1:end,np+1:end) ;
int_function_matrix = int_function(p,t) ;

% linear solution with Ru = Rv/rq = Vmu/Kmu*Cu
Cu_init = (Ku + (variables.Vmu/variables.Kmu).*int_function_matrix - Kbu)\(Bu);
Cv_init = (Kv - Kbv)\(variables.rq.*(variables.Vmu/variables.Kmu).*int_function_matrix*Cu_init+Bv);

% initialize solution for newton-raphson
Cu = Cu_init ; Cv = Cv_init ;

%NEWTON-RAPHSON EXPLICIT SOVER
h = waitbar(0,'Initializing') ; % creating a waitbar

for it = 1:max_it
    % Compute source term and its jacobian
    [Ru,Rv,RudCu,RudCv,RvdCu,RvdCv] = fun_diff([Cu; Cv], variables) ;
    
    % System function
    F = [Ku*Cu + int_function_matrix*Ru - Kbu*Cu - Bu ; ...
        Kv*Cv - int_function_matrix*Rv - Kbv*Cv - Bv] ;
    
    % Jacobian of system
    JF = [[Ku + int_function_matrix*RudCu - Kbu,  int_function_matrix*RudCv]; ...
        [-int_function_matrix*RvdCu,             Kv - int_function_matrix*RvdCv - Kbv]];
    
    sol = [Cu; Cv] - JF\F ;
    
    % verify tolerance
    tol = norm(sol-[Cu; Cv]) ;
    if tol < tol_min ; break ; end
    
    % update buffer
    Cu = sol(1:np) ;
    Cv = sol(np+1:end) ;
    
    % update waitbar
    if mod(it,10)==0
        waitbar(it/max_it,h,['tol = ' num2str(tol) newline ...
            'iteration = ' num2str(it) ' / ' num2str(max_it)]) ;
    end
end

delete(h) ;

%% VIEW SOL
if plot_sol && strcmp(type,'pear')
    figure ; hold on ;
    
    [xq,yq] = meshgrid(linspace(min(p(:,1)),max(p(:,1)),200), ...
        linspace(min(p(:,2)),max(p(:,2)),200));
    
    subplot(1,2,1) ;
    title('O_2 concentration') ;
    gsolu = griddata(p(:,1),p(:,2),full(Cu),xq,yq) ;
    [~,h] = contourf(xq,yq,gsolu) ;
    set(h,'LineColor','none') ;
    % colorbar ;
    xlabel('O_2') ;
    colormap('parula') ;
    caxis([0 10]) ;
    set(gca,'xtick',[]) ;
    set(gca,'ytick',[])
    
    subplot(1,2,2) ;
    title('CO_2 concentration') ;
    gsolv = griddata(p(:,1),p(:,2),full(Cv),xq,yq) ;
    [~,h] = contourf(xq,yq,gsolv) ;
    set(h,'LineColor','none') ;
    colorbar ;
    xlabel('CO_2') ;
    colormap('parula') ;
    caxis([0 10]) ;
    set(gca,'xtick',[]) ;
    set(gca,'ytick',[]) ;
    
    suptitle([casus newline ...
        'T = ' num2str(Tcel) '°C / ' ...
        'O_2 = ' num2str(eta_u) '% / ' ...
        'CO_2 = ' num2str(eta_v) '%']) ;  
end

%% RETURN (not ready for the moment)
assert(nargout==2, 'Wrong number of output arguments (2)') ;
varargout{1} = Cu ;
varargout{2} = Cv ;

if strcmp(type,'test')
    R = max(p(:,1)) ; % radius of test mesh
    
    r_middle_idx = and(p(:,2)<1e-3, p(:,2)>-1e-3) ;
    r_middle = p(r_middle_idx,1) ;
    [r_sort, idx_sort] = sort(r_middle) ;
    
    param = [Tcel ; eta_u ; eta_v] ;
    odefun_param = @(x,y) odefun(x,y,param) ;
    bcfun_param = @(x,y) bcfun(x,y,param) ;
    
    solinit = bvpinit(linspace(0,R,50), [0; 1]) ;
    bvp_opt = bvpset('SingularTerm',[0 0 ; 0 -1]) ;
    Cu_init_real = bvp5c(odefun_param, bcfun_param, solinit, bvp_opt) ;
    Cu_init_real = interp1(Cu_init_real.x,Cu_init_real.y(1,:),r_sort') ;
    
    Cu_init = Cu_init(r_middle_idx) ;
    Cu_init = Cu_init(idx_sort) ;
    
    figure ; hold on ;
    plot(r_sort,full(Cu_init(idx_sort)),'-k','LineWidth',2) ;
    plot(r_sort,Cu_init_real',':k','LineWidth',2) ;
    title([casus newline 'Approximation error']) ;
    legend({'FEM model', 'Real solution'}) ;
    xlabel('r') ; ylabel('C_u') ;
    %axis([0 R, 0 1]) ;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dydx = odefun(~,y,param)

Tcel = param(1) ;
eta_u = param(2) ;
eta_v = param(3) ;

variables = generate_variables('pear',Tcel, eta_u, eta_v) ;

% dydx = [y(2) ; -y(2)/x + (variables.Vmu/variables.Kmu)*y(1)/variables.Dur] ;
dydx = [y(2) ; (variables.Vmu/variables.Kmu)*y(1)/variables.Dur] ;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res = bcfun(ya,yb,param)

Tcel = param(1) ;
eta_u = param(2) ;
eta_v = param(3) ;

variables = generate_variables('pear',Tcel, eta_u, eta_v) ;

res = [ya(2) ; variables.hu*(yb(1)-variables.Cuamb)/variables.Dur+yb(2)] ;

end