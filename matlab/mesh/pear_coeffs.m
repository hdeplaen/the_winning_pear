function [P_tot] = pear_coeffs()
%PEAR_COEFFS returns the bezier coefficients of the pear
%   Author: Henri De Plaen

scale_x = .3 ;
scale_y = 1 ;

P1 = [305,125.8 ; ...
    332.7,125.8 ; ...
    346.6,164.9 ; ...
    346.4,208.8] ;

P2 = [346.4,208.8 ; ...
    346.2,241.7 ; ...
    339.7,246.7 ; ...
    356.8, 279.7] ;

P3 = [356.8, 279.7 ; ...
    369.4, 304 ; ...
    386.8, 301.1 ; ...
    395, 336.5 ] ;

P4 = [395, 336.5 ; ...
    402.2, 367.7 ; ...
    395.6, 382.9 ; ...
    374.1, 405.6 ] ;

P5 = [374.1, 405.6 ; ...
    345.2, 436.1 ; ...
    336.7, 435.7 ; ...
    305, 436.5 ] ;

P6 = [305, 436.5 ; ...
    305, 436.5 ; ...
    305,125.8 ; ...
    305,125.8 ] ;    

P = [P1 ; P2 ; P3 ; P4 ; P5] ;
max_P_x = max(max(P(:,1))) ; min_P_x = min(min(P(:,1))) ;
max_P_y = max(max(P(:,2))) ; min_P_y = min(min(P(:,2))) ;
diff_P_x = (max_P_x-min_P_x)/scale_x ;
diff_P_y = (max_P_y-min_P_y)/scale_y ;

P1(:,1) = (P1(:,1)-min_P_x)./diff_P_x ;
P2(:,1) = (P2(:,1)-min_P_x)./diff_P_x ;
P3(:,1) = (P3(:,1)-min_P_x)./diff_P_x ;
P4(:,1) = (P4(:,1)-min_P_x)./diff_P_x ;
P5(:,1) = (P5(:,1)-min_P_x)./diff_P_x ;
P6(:,1) = (P6(:,1)-min_P_x)./diff_P_x ;

P1(:,2) = 1-(P1(:,2)-min_P_y)./diff_P_y ;
P2(:,2) = 1-(P2(:,2)-min_P_y)./diff_P_y ;
P3(:,2) = 1-(P3(:,2)-min_P_y)./diff_P_y ;
P4(:,2) = 1-(P4(:,2)-min_P_y)./diff_P_y ;
P5(:,2) = 1-(P5(:,2)-min_P_y)./diff_P_y ;
P6(:,2) = 1-(P6(:,2)-min_P_y)./diff_P_y ;

P_tot(:,:,1) = P1 ;
P_tot(:,:,2) = P2 ;
P_tot(:,:,3) = P3 ;
P_tot(:,:,4) = P4 ;
P_tot(:,:,5) = P5 ;
P_tot(:,:,6) = P6 ;

end