function variables = generate_variables(type,Tcel, eta_u, eta_v)
%Generates mesh for the pear problem
%   variables = generate_mesh(type)
%
% INPUT
%   type   : 'pear' or 'test'
%   Tcel   : parameter 1
%   eta_u  : parameter 2
%   eta_v  : parameter 3
% OUTPUT
%   variables: problem parameters
%
%Author: Henri De Plaen
%KU Leuven
%Project WIT: the winning pear
%Date: Nov 2018

assert(nargin==4, 'Wrong number of input arguments (4)') ;
assert(nargout==1, 'Wrong number of output arguments (1)') ;

switch type
    case 'pear' ; variables = var_pear(Tcel, eta_u, eta_v) ; 
    case 'test' ; variables = var_test() ;
    otherwise ; error('Type not recongized') ;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = var_pear(Tcel, eta_u, eta_v) 

assert(nargin==3, 'No input arguments needed') ;

variables = struct ;

%% PARAMS
variables.Dur = 2.8e-10 ;
variables.Duz = 1.1e-9 ;
variables.Dvr = 2.32e-9 ;
variables.Dvz = 6.97e-9 ;

variables.Tcel = Tcel ;
variables.eta_u = eta_u/100 ;
variables.eta_v = eta_v/100 ;

variables.Rg = 8.314 ;
variables.hu = 7e-7 ;
variables.hv = 7.5e-7 ;
variables.patm = 101300 ;
variables.T = variables.Tcel+273.15 ;
variables.Tref = 20 + 273.15 ;

variables.Cuamb = variables.patm*variables.eta_u/variables.Rg/variables.T ;
variables.Cvamb = variables.patm*variables.eta_v/variables.Rg/variables.T ;

variables.rq = 0.97 ;
variables.Kmfu = 0.1149 ;
variables.Kmv = 27.2438 ;
variables.Kmu = 0.4103 ;

variables.Eavmfvref = 56700 ;
variables.Vmfvref = 1.61e-4 ;
variables.Vmfv = variables.Vmfvref*exp(variables.Eavmfvref/variables.Rg* ...
    (1/variables.Tref-1/variables.T)) ;

variables.Eavmuref = 80200 ;
variables.Vmuref = 2.39e-4 ;
variables.Vmu = variables.Vmuref*exp(variables.Eavmuref/variables.Rg* ...
    (1/variables.Tref-1/variables.T)) ;

%% RETURN
assert(nargout==1, 'Outuput is only one structure argument') ;
varargout{1} = variables ;

end

function varargout = var_test() 

assert(nargin==0, 'Wrong number of input arguments (3)') ;

variables = struct ;

%% PARAMS
variables.Dur = 1 ;
variables.Duz = 1 ;
variables.Dvr = 1 ;
variables.Dvz = 1 ;

variables.Tcel = 1 ;
variables.eta_u = 1 ;
variables.eta_v = 1 ;

variables.Rg = 1 ;
variables.hu = 1 ;
variables.hv = 1 ;
variables.patm = 1 ;
variables.T = 1 ;
variables.Tref = 1 ;

variables.Cuamb = 1 ;
variables.Cvamb = 1 ;

variables.rq = 1 ;
variables.Kmfu = 1 ;
variables.Kmv = 1 ;
variables.Kmu = 1 ;

variables.Eavmfvref = 1 ;
variables.Vmfvref = 1 ;
variables.Vmfv = 1 ;

variables.Eavmuref = 1 ;
variables.Vmuref = 1 ;
variables.Vmu =  1 ;

%% RETURN
assert(nargout==1, 'Outuput is only one structure argument') ;
varargout{1} = variables ;

end