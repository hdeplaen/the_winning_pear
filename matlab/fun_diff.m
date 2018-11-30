function [Ru,Rv,RudCu,RudCv,RvdCu,RvdCv] = fun_diff(sol, variables)
%Evalute the source term and its jacobian at given points
%   [Ru,Rv,RudCu,RudCv,RvdCu,RvdCv] = fun_diff(sol, variables)
%
% INPUT
%   sol: nodes solution at which evaluate the function
%   variables: problem parameters
% OUTPUT
%   source term and its jacobian
%
%Author: Henri De Plaen
%KU Leuven
%Project WIT: the winning pear
%Date: Nov 2018

assert(nargin==2, 'Wrong number of input arguments (2)') ;

%% PRELIMINARIES
n_sol = size(sol,1) ;           % number of points (2*np normally)

assert(mod(n_sol,2)==0) ;       % verifying that the number of points if even
Cu = sol(1:n_sol/2,:) ;         % sub-vector corresponding to Cu
Cv = sol(n_sol/2+1:end,:) ;     % sub-vector corresponding to Cv

%% COMPUTE SOURCE TERM
% SOURCE TERM FUNCTION
Ru = -variables.Vmu.*Cu./((variables.Kmu+Cu).*(1+Cv./variables.Kmv)) ;

Rv = -variables.rq.*Ru + variables.Vmfv./(1+Cu./variables.Kmfu) ;

% SOURCE TERM JACOBIAN
RudCu = -(variables.Kmu.*variables.Kmv.*variables.Vmu)./ ...
    ((Cu + variables.Kmu).^2.*(Cv + variables.Kmv)) ;

RudCv = (Cu.*variables.Kmv.*variables.Vmu)./ ...
    ((Cu + variables.Kmu).*(Cv + variables.Kmv).^2) ; 

RvdCu = (variables.Kmv.*variables.Vmu.*variables.rq)./ ...
    ((Cu + variables.Kmu).*(Cv + variables.Kmv)) - ...
    (variables.Kmfu.*variables.Vmfv)./(Cu + variables.Kmfu).^2 - ...
    (Cu.*variables.Kmv.*variables.Vmu.*variables.rq)./ ...
    ((Cu + variables.Kmu).^2.*(Cv + variables.Kmv)) ;

RvdCv = -(Cu.*variables.Kmv.*variables.Vmu.*variables.rq)./ ...
    ((Cu + variables.Kmu).*(Cv + variables.Kmv).^2) ;

end
