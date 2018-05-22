
function cov = build_cov( ptmp,coord_float,po_system_configuration);

% This function was created from build_ptmp_cov  which the description is written below:
% 
% Description of build_ptmp_cov 
% ____
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
%____
%
% Cecile Cabanes, June 2013 
% Compared to build_ptmp_cov, build_cov( ptmp,coord_float,po_system_configuration) take into account the horizontal covariance between different mapped profiles. 
% This lateral covariance take into account the fact that a mapped proﬁle on a Argo position is build from a set of historical proﬁles that is not very diﬀerent from the set used to build a mapped proﬁle at the next or previous Argo proﬁle position.
% This lateral covariance between two mapped proﬁles is constructed using a Gaussian function and the large spatial scales
%
% Cecile Cabanes, 2017
% use the small spatial scales instead of the large spatial scales to build
% the lateral covariance: found to be the best compromise between
% informative errors  and large enough NDF for AIC criterium, at least for
% the scales defined for the North Atlantic bassin

% Set up the theta boundaries for water masses

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

cov=zeros(m*n,m); % set up the output covariance matrix, with 1's along the diagonal


for k=1:n % for each profile
    k1 = (k-1)*m;
    for l=1:1 % for each profile
        l1 = (l-1)*m;
        for i=1:m % for all levels
            for j=1:m
                if(i<j) % upper triangle, look down the water column for vertical scale
                    Ltheta = interp1( ptboundaries, ptscale_down, ptmp(i,k), 'linear' );
                    cov(i+k1,j+l1) = exp( - ( ptmp(j,k) - ptmp(i,k) ).^2/ Ltheta.^2 );
                elseif(i>j) % lower triangle, look up the water column for vertical scale
                    Ltheta = interp1( ptboundaries, ptscale_up, ptmp(i,k), 'linear' );
                    cov(i+k1,j+l1) = exp( - ( ptmp(j,k) - ptmp(i,k) ).^2/ Ltheta.^2 );
                elseif(i==j)   
                     cov(i+k1,j+l1) =1;
                end
                if isnan(cov(i+k1,j+l1))
                    cov(i+k1,j+l1)=1;
                end
            end
        end
    end
end


% 
% set up the horizontal covariance matrix  
scale_lon = str2num(po_system_configuration.MAPSCALE_LONGITUDE_SMALL);
scale_lat = str2num(po_system_configuration.MAPSCALE_LATITUDE_SMALL);
scale_phi = str2num(po_system_configuration.MAPSCALE_PHI_SMALL);
map_use_pv = str2num(po_system_configuration.MAP_USE_PV);

for k=1:n
  b(k,:)=covarxy_pv(coord_float(k,:),coord_float,scale_lon,scale_lat,scale_phi,map_use_pv);
end
b=b(:,1:n);

covh=b;

%keyboard
% build the final covariance matrix, taking into account horizontal and vertical covariance

covn=repmat(cov,[1,n]);
for k=1:n % for each profile
    kd = (k-1)*m+1;
    kf =  k*m;
    for l=1:n % for each profile
        ld = (l-1)*m+1;
        lf = (l)*m;
        covn(kd:kf,ld:lf)=covh(k,l)*covn(kd:kf,ld:lf);
    end
end

cov=covn;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [a] = covarxy_pv(x1,x2,gs,ts,phi,pv)

function  [a] = covarxy_pv(x1,x2,gs,ts, phi, pv)

m=size(x1);
n=size(x2);
a=zeros(m(1),n(1));
one=ones(n(1),1);


% derive planetary vorticity at each data point ----

Zx1=x1(:,3); % depth at x1 datapoint
Zx2=x2(:,3); % depth at x2 datapoint

PV_x1=(2*7.292*10^-5.*sin(x1(:,1).*pi/180))./Zx1;
PV_x2=(2*7.292*10^-5.*sin(x2(:,1).*pi/180))./Zx2;


% calculate covariance term, use pv as optional ----

for i=1:m(1)
   tmp=((x1(i,1)-x2(:,1))./ts).^2 + ...
       ((x1(i,2)-x2(:,2))./gs).^2;
   if(pv &PV_x1~=0&PV_x2~=0)  % if PV is wanted and the point is not on the equator ---
       tmp = tmp + ...
             ( (PV_x1(i)-PV_x2)./sqrt( PV_x1(i).^2+PV_x2.^2 )./phi ).^2;
   end

  
   a(i,:)=exp(-tmp');
end

return
