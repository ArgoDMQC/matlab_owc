
function [noisevariance]=noise(a,xlat,xlong)

% function [noisevariance]=noise(a,xlat,xlong)
%
% INPUT:
% xlat = a column vector (nx1)
% xlong = a column vector (nx1)
% a = m*n matrix fields consisting of m layers and n casts
% - this can contain NaN's.
%
% OUTPUT:
% noisevariance = variance of noise between closest casts
% - an estimate of the noise due to variations between casts,
%   one variance value per level.
%
% NOTE:
% noise = measured value - true value.
% 2*noisevariance = variance between casts. Here, the variance is
% based on closest casts.

[m,n]=size(a);
diff=NaN*ones(m,n);
noisevariance=zeros(m,1);
for i=1:n
	xlat0=xlat(i)*ones(n,1);
	xlong0=xlong(i)*ones(n,1);
	r=(xlat0-xlat).^2+ (xlong0-xlong).^2;
	index=find( r > 0);
	[tmp]=min(r(index));
	j=find(r==tmp);
	diff(:,i)=a(:,i)-a(:,j(1));
end
for i=1:m
	ii=find(isnan(diff(i,:))==0);
	noisevariance(i)=sum(diff(i,ii).^2)/(2*length(ii));
end

