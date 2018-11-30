function variables = generate_variables(type)
%Generates mesh for the pear problem
%   variables = generate_mesh(type)
%
% INPUT
%   type: 'pear' or 'test'
% OUTPUT
%   variables: problem parameters
%
%Author: Henri De Plaen
%KU Leuven
%Project WIT: the winning pear
%Date: Nov 2018

assert(nargin==1, 'Wrong number of input arguments (1)') ;
assert(nargout==1, 'Wrong number of output arguments (1)') ;

switch type
    case 'pear' ; variables = var_pear() ; 
    case 'test' ; variables = var_test() ;
    otherwise ; error('Type not recongized') ;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = var_pear() 

assert(nargin==0, 'No input arguments needed') ;

variables = struct ;

%% PARAMS
variables.Dur = 2.8e-10 ;
variables.Duz = 1.1e-9 ;
variables.Dvr = 2.32e-9 ;
variables.Dvz = 6.97e-9 ;

variables.Tcel = 25 ;
variables.eta_u = 20.8/100 ;
variables.eta_v = 0.04/100 ;

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

assert(nargin==0, 'No input arguments needed') ;

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
% variables.Cuamb = variables.patm*variables.eta_u/variables.Rg/variables.T ;
% variables.Cvamb = variables.patm*variables.eta_v/variables.Rg/variables.T ;

variables.rq = 1 ;
variables.Kmfu = 1 ;
variables.Kmv = 1 ;
variables.Kmu = 1 ;

variables.Eavmfvref = 1 ;
variables.Vmfvref = 1 ;
variables.Vmfv = 1 ;
% variables.Vmfv = variables.Vmfvref*exp(variables.Eavmfvref/variables.Rg* ...
%     (1/variables.Tref-1/variables.T)) ;

variables.Eavmuref = 1 ;
variables.Vmuref = 1 ;
variables.Vmu =  1 ;
% variables.Vmu = variables.Vmuref*exp(variables.Eavmuref/variables.Rg* ...
%     (1/variables.Tref-1/variables.T)) ;

%% RETURN
assert(nargout==1, 'Outuput is only one structure argument') ;
varargout{1} = variables ;

end