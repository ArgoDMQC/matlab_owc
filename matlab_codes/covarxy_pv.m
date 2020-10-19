
% function [a] = covarxy_pv(x1,x2,gs,ts,phi,pv)

% This function is based on covarxyt_pv.m
% Output:
% a 	m*n  unscaled Gaussian covariance function
%
% Input:
% x1	m*4 model grid, col 1+2+3+4 are lat, long, date, depth respectively
% x2	n*4 data grid, col 1+2+3+4 are lat, long, date, depth respectively
% ts	scalar, horizontal latitude scale  (degrees)
% gs    scalar, horizontal longitude scale (degrees)
% phi   scalar cross isobaric scale (for depth dependence)
% pv    include potential vorticity in the covariance or not (1=include, 0=exclude)
%
% Breck Owens, modified 21 November 2007
% A. Wong, 29 April 2005

function  [a] = covarxy_pv(x1,x2,gs,ts, phi, pv)

m=size(x1);
n=size(x2);
a=zeros(m(1),n(1));
one=ones(n(1),1);


% derive planetary vorticity at each data point ----

Zx1=x1(:,4); % depth at x1 datapoint
Zx2=x2(:,4); % depth at x2 datapoint

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



