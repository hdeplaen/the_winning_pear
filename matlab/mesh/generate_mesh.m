% Author: Henri De Plaen
% KULeuven
% Project WIT : winning pear
% Date: Nov 2018

clear all ; close all ; clc ;

%% PDE geom and mesh
model = createpde(1) ;
geom = geometryFromEdges(model,@pear_geometry) ;
mesh = generateMesh(model,'GeometricOrder','linear', 'Hgrad', 1.9, 'Hmax', 0.2);
[p,e,t] = meshToPet(mesh);
outerpoints = length(e) ;
innerpoints = length(p) - outerpoints ;

%% plot
pdemesh(model,'NodeLabels','on') ; %,'ElementLabels','on') ;
title('Pear mesh') ;
xlabel('x') ; ylabel('y') ;
axis([-.1 0.35 -.1 1.1]) ;

%% export csv
node_csv = t(1:3,:)'-1 ;
xp_csv = p(1,:)' ;
yp_csv = p(2,:)' ;

csvwrite('node.csv',node_csv) ;
csvwrite('xp.csv',xp_csv) ;
csvwrite('yp.csv',yp_csv) ;


p = p' ;
e = e' ;
t = t' ;