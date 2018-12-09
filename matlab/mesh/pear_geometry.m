function [ x,y ] = pear_geometry( bs,s )
%PEAR Creates the geometry function of a pear
%   Author: Henri De Plaen
P = pear_coeffs() ;

% dl = [P(1,2,1),P(1,2,2),P(1,2,3),P(1,2,4),P(1,2,5),P(1,2,6); % start parameter values
%     P(4,2,1),P(4,2,2),P(4,2,3),P(4,2,4),P(4,2,5),P(4,2,6); % end parameter values
%     1,1,1,1,1,1; % region label to left
%     0,0,0,0,0,0]; % region label to right
    
dl = [0,0,0,0,0,0 ; 1,1,1,1,1,1 ; 0,0,0,0,0,0 ; 1,1,1,1,1,1 ] ;

switch nargin
    case 0
        x = 6; % six edge segments
        return
        
    case 1
        x = dl(:,bs); % return requested columns
        return
        
    case 2
        [m,n]=size(bs);
        if m==1 && n==1
            bs=bs*ones(size(s)); % expand bs
        elseif m~=size(s,1) || n~=size(s,2)
            error('bs must be scalar or of same size as s');
        end
        
        npoints = 400 ;
        
        x = [] ;
        y = [] ;
        
        bs_vec = bs(:) ;
        s_vec = s(:) ;
        for nbs = 1:max(bs_vec)
            idx = bs_vec==nbs ;
            s_loc = s_vec(idx) ;
            [x_loc, y_loc] = BezierCurve(linspace(0,1,npoints),P(:,:,nbs)) ;
            s_loc=pdearcl(linspace(0,1,npoints),[x_loc(:)';y_loc(:)'],s_loc(:),0,1);
            [x_loc, y_loc] = BezierCurve(s_loc',P(:,:,nbs)) ;
            x = [x ; x_loc] ;
            y = [y ; y_loc] ;
        end
        
        x = reshape(x,size(bs)) ;
        y = reshape(y,size(bs)) ;
end
end

