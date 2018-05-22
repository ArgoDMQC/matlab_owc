
function pa_decimal_dates = cal2dec(pa_month, pa_day, pa_hour, pa_minute)

%function t= cal2dec(pa_month, pa_day, pa_hour,pa_minute)
%
%converts calander pa_day to decimal pa_day
%
%P. Robbins 94

if nargin < 3
  pa_hour = 0 ;
end
if nargin < 4
  pa_minute = 0 ;
end

la_months = [0 31 28 31 30 31 30 31 31 30 31 30 31];
ln_cumulative_months = cumsum(la_months);


pa_decimal_dates = ln_cumulative_months( pa_month )' + pa_day - 1 + pa_hour/24 + pa_minute/24/60;

