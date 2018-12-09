%% PRELIMIARIES
delete(get(0,'Children')) ;     % delete restant waitbars
clear all ; close all ; clc ;   % start again from scratch

% %% RUN TESTS
% hmax = 0.0008 ;
% type = 'test' ;
% refer = 'radial' ;
% 
% [~, ~] = pear_problem(hmax,type,1,refer,-1,2,0.7,'Optimal CA') ;
% [~, ~] = pear_problem(hmax,type,1,refer,-1,2,5,'Disorder inducing') ;
% [~, ~] = pear_problem(hmax,type,1,refer,-1,20.8,0,'Pre-cooling') ;
% [~, ~] = pear_problem(hmax,type,1,refer,7,20.8,0,'Refrigerator') ;
% [~, ~] = pear_problem(hmax,type,1,refer,20,20.8,0,'Shelf life') ;

%% RUN PEARS
hmax = 0.0008 ;
type = 'pear' ;
refer = 'radial' ;
% 
% [~, ~] = pear_problem(hmax,type,1,refer,-1,2,0.7,'Optimal CA') ;
% [~, ~] = pear_problem(hmax,type,1,refer,-1,2,5,'Disorder inducing') ;
% [~, ~] = pear_problem(hmax,type,1,refer,-1,20.8,0,'Pre-cooling') ;
% [~, ~] = pear_problem(hmax,type,1,refer,7,20.8,0,'Refrigerator') ;
[~, ~] = pear_problem(hmax,type,1,refer,20,20.8,0,'Shelf life') ;