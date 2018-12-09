function ep = eval_geometry(geom_function, prec)

ne = feval(geom_function) ;                     % number of edges

ep = zeros(ne*prec,2) ;                      % prealloc
s = linspace(0,1,prec) ;                       % arc range

for idx = 1:ne
    [x,y] = feval(geom_function,idx,s) ;
    ep( (idx-1)*prec+1 : (idx)*prec, 1) = x(:) ;
    ep( (idx-1)*prec+1 : (idx)*prec, 2) = y(:) ;
end

end