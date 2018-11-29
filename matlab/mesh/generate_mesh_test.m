% Author: Henri De Plaen
% KULeuven
% Project WIT : winning pear
% Date: Nov 2018

clear all ; close all ; clc ;

%% PDE geom and mesh
model = createpde(1) ;

C1 = [1,0,1,1]';
R1 = [3,4,-1,-1,0,0,2,0,0,2]' ;
C1 = [C1;zeros(length(R1) - length(C1),1)];
gm = [R1,C1];
sf = 'C1-R1';
ns = char('R1','C1');
ns = ns';
g = decsg(gm,sf,ns);

geom = geometryFromEdges(model,g) ;
mesh = generateMesh(model,'GeometricOrder','linear', 'Hgrad', 1.9, 'Hmax', 0.15);
[p,e,t] = meshToPet(mesh);
outerpoints = length(e) ;
innerpoints = length(p) - outerpoints ;

%% plot
pdemesh(model,'NodeLabels','on') ; %,'ElementLabels','on') ;
title('Test mesh') ;
xlabel('x') ; ylabel('y') ;
axis([-.1 1.1 -.1 2.1]) ;

%% export csv
node_csv = t(1:3,:)'-1 ;
xp_csv = p(1,:)' ;
yp_csv = p(2,:)' ;

csvwrite('node_test.csv',node_csv) ;
csvwrite('xp_test.csv',xp_csv) ;
csvwrite('yp_test.csv',yp_csv) ;


p = p' ;
e = e' ;
t = t' ;