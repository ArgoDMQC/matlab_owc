
function cov = build_ptmp_cov( ptmp );

%
% This function builds the covariance matrix (a square matrix) that has nxn number
% of tiles, each tile is of mxm size.
%
% The vertical covariance matrix is the building tile. It has 1 at its diagonal,
% then decreases exponentially in the off-diagonals, which represents the vertical
% covariance between water masses.
%
% This version assums that each profile is independent from the other ones.
%
% Annie Wong, May 2001
% Breck Owens, November 2007
%

% set up the theta boundaries for water masses

ptboundaries = [ 30; 24; 18; 12; 8; 4; 2.5; 1; -2 ];
ptscale_down = [ 6; 6; 6; 4; 4; 1.5; 1.5; 1; 1 ];
ptscale_up = [ 6; 6; 6; 6; 4; 4; 1.5; 1.5; 1 ];


% set up the building tile = vertical covariance matrix
%
% upper triangle of the matrix = covariance of each ptlevel with every ptlevel below it,
%  looking down the water column from the diagonal
% lower triangle of the matrix = covariance of each ptlevel with every ptlevel above it,
%  looking up the water column from the diagonal

[m,n]=size(ptmp);

cov=eye(m*n,m*n); % set up the output covariance matrix, with 1's along the diagonal

for k=1:n % for each profile
    k1 = (k-1)*m;
    for i=1:m % for all levels
        for j=1:m
            if(i<j) % upper triangle, look down the water column for vertical scale
                Ltheta = interp1( ptboundaries, ptscale_down, ptmp(i,k), 'linear' );
                cov(i+k1,j+k1) = exp( - ( ptmp(j,k) - ptmp(i,k) ).^2/ Ltheta.^2 );
            elseif(i>j) % lower triangle, look up the water column for vertical scale
                Ltheta = interp1( ptboundaries, ptscale_up, ptmp(i,k), 'linear' );
                cov(i+k1,j+k1) = exp( - ( ptmp(j,k) - ptmp(i,k) ).^2/ Ltheta.^2 );
            end
        end
    end
end

return
