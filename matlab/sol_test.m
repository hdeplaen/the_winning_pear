% Author: Henri De Plaen

clear all ; close all ; clc ;

generate_mesh ; close all ;
np = size(p,1) ;

regv_init = 1e+2 ;
regv = 1 ;

h = waitbar(0,'Initializing') ;
max_it = 1e+5 ;
tol_min = 1e-10 ;

K = stiff_matrix(p,t) ;
e6 = e(e(:,5)==6,:) ;
ibcd = union(e(:,1),e(:,2)) ;
ibcd6 = union(e6(:,1),e6(:,2)) ;
%K([ibcd6 ibcd6+np],[ibcd6 ibcd6+np]) = K([ibcd6 ibcd6+np],[ibcd6 ibcd6+np]) + 10^15*speye(2*length(ibcd6)) ;

[B,Rb] = boundary_vector_bis(p,e) ;

np = size(p,1) ;
sol_buff = abs(.01*randn(2*np,1)) ;
sol = zeros(2*np,1) ; %prealloc
tol = 1 ;

[F,dF] = fun_vector(p,t,sol_buff) ;
F([ibcd ibcd+np]) = 0 ;
dF([ibcd ibcd+np],[ibcd ibcd+np]) = 0 ;
dF = 1/regv_init*dF ;
sol = (K-Rb-dF)\(F-dF*sol_buff+B) ;

for it = 1:max_it
    [F,dF] = fun_vector(p,t,sol_buff) ;
    F([ibcd ibcd+np]) = 0 ;
    dF([ibcd ibcd+np],[ibcd ibcd+np]) = 0 ;
    dF = 1/regv*dF ;
    
    %     K = full(K) ;
    %     F = full(F) ;
    %     B = full(B) ;
    %     Rb = full(Rb) ;
    %     dF = full(dF) ;
    
    sol = (K-Rb-dF)\(F-dF*sol_buff+B) ;
    
    tol = sqrt(mean((sol-sol_buff).^2)) ;
    if abs(tol) < tol_min
        break ;
    end
    sol_buff = sol ;
    
    if mod(it,300)==0
        waitbar(it/max_it,h,['tol = ' num2str(tol)]) ;
    end
end

delete(h) ;

%% VIEW SOL
sol_u = sol(1:np) ;
sol_v = sol(np+1:end) ;

figure ; hold on ;
dm = diag(K) ;

subplot(1,2,1) ;
title('Oxygen concentration') ;
pdeplot(model,'XYData',sol_u(:)) ;
%pdeplot(model,'XYData',dm(1:np)) ;
xlabel('r') ; ylabel('z') ;
colormap('parula') ;

subplot(1,2,2) ;
title('Carbon dioxide concentration') ;
pdeplot(model,'XYData',sol_v(:)) ;
%pdeplot(model,'XYData',dm(np+1:end)) ;
xlabel('r') ; ylabel('z') ;
colormap('parula') ;