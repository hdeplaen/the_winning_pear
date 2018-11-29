function varargout = var_test() 

assert(nargin==0, 'No input arguments needed') ;

variables = struct ;
%% PARAMS (unused variables in comment)
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