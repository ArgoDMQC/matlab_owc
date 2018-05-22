
function [signalvariance]=signal(a)

% function [signalvariance]=signal(a)
%
% INPUT:
% a = m*n matrix (m layers, n casts) field.
% - this can contain NaN's
%
% OUTPUT:
% signalvariance = an estimate of variance of signal at each depth level.
%
% NOTE:
% signal = data - bestfit(data). Here, bestfit is simply the average.

[m,n]=size(a);
signal=NaN*ones(m,n);
signalvariance=zeros(m,1);
for i=1:m
	jj=find(isnan(a(i,:))==0);
	signal(i,jj)=a(i,jj)-mean(a(i,jj));
	signalvariance(i)=(std(signal(i,jj))^2)*(length(jj)-1)/length(jj);
end

% signalvariance(i)=(signal(i,jj)*signal(i,jj)')/length(jj)
% will give the same results, because mean(signal)=0.


