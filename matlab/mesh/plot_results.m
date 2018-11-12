% Author: Henri De Plaen
% KULeuven
% Project WIT : pear
% Date: March 2018

generate_mesh ;
close all ; clc ;

%% LOAD RESULTS
sol = csvread('../src/sol.csv') ;
lsol = length(sol)/2 ;
sol_u = sol(1:lsol) ;
sol_v = sol(lsol+1:end) ;

%% PLOT
figure ; hold on ;

subplot(1,2,1) ;
title('Oxygen concentration') ;
pdeplot(model,'XYData',sol_u(:)) ;
xlabel('r') ; ylabel('z') ;

subplot(1,2,2) ;
title('Carbon dioxide concentration') ;
pdeplot(model,'XYData',sol_v(:)) ;
xlabel('r') ; ylabel('z') ;

%% MISC
stiff = csvread('../src/stiff.csv') ;