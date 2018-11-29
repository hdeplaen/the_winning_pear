function varargout = var_real() 

assert(nargin==0, 'No input arguments needed') ;

variables = struct ;
%% PARAMS (unused variables in comment)
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