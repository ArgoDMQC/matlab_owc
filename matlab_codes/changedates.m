
% This function changes dates in format YYYYMMDDhhmmss to a
% decimal number. The input dates can be in either a single
% row or a single column, but the output dates are organised
% in a single column.
%
% A. Wong, 29 May 2001
%

function [dates]=changedates(dates_format2);



dates=NaN.*ones(length(dates_format2),1); %organise dates in a single column

for i=1:length(dates_format2)
 if(isnan(dates_format2(i))==0)
  junk=int2str(dates_format2(i));
  yr=str2num(junk(:,1:4));
  mo=str2num(junk(:,5:6));
  day=str2num(junk(:,7:8));
  hr=str2num(junk(:,9:10));
  min=str2num(junk(:,11:12));
  if(mo<1|mo>12|day<1|day>31)
   dates(i)=yr;
  else
   dates(i)=yr+cal2dec(mo,day,hr,min)./365;
  end
 end
end

