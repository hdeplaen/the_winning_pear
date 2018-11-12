function [x, y] = BezierCurve(u, P) 
% This function constructs a Bezier curve from given control points. P is a 
% vector of control points. u is the vector of points to calculate. 
% 
% Example: 
% 
% P = [0 0; 1 1; 2 5; 5 -1]; 
% y = BezierCurve(linspace(0,5,1000), P); 
% plot(x, y, P(:, 1), P(:, 2), 'x-', 'LineWidth', 2); set(gca, 'FontSize', 16) 
% 
% Prakash Manandhar, pmanandhar@umassd.edu
% modif: Henri De Plaen

u = u(:) ;
Np = size(P, 1);  
B = zeros(length(u), Np); 
for i = 1:Np 
B(:,i) = nchoosek(Np-1,i-1).*u.^(i-1).*(1-u).^(Np-i); %B is the Bernstein polynomial value 
end 
S = B*P; 
x = S(:, 1); 
y = S(:, 2);

x = x(:) ;
y = y(:) ;

end