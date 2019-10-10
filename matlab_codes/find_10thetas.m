
function [tlevels, plevels, index, var_s_Thetalevels, Thetalevels] = find_10thetas( SAL, PTMP, PRES, la_ptmp, use_theta_gt, use_theta_lt, use_pres_gt, use_pres_lt, use_percent_gt )

%-----------------------------------------------------------------------------------------
% Chooses 10 theta levels from the float series for use in the linear fit.
% These 10 theta levels are the ones with the minimum S variance on theta.
%
% These 10 theta levels are distinct (ie. they don't repeat each other).
%
% OUTPUT:
% tlevels = the 10 chosen theta levels
% plevels = nominal pressure of the 10 chosen theta levels
% index = indices corresponding to the float_source matrices of the 10 chosen theta levels
% var_s_Thetalevels = salinity variance on thetas
% Thetalevels = thetas on which the salinity variances are calculated
%
% A.Wong, March 2009.
% Giulio Notarstefano (GN), Dec 2018: allows theta min and max values with a 2-digits precision
%------------------------------------------------------------------------------------------

%SAL=unique_SAL;
%PTMP=unique_PTMP;
%PRES=unique_PRES;
%la_ptmp=unique_la_ptmp;


% initialize variables ---

h=10; % chooses 10 theta levels
[m,n]=size(PRES);
index=NaN.*ones(h,n);
tlevels=NaN.*ones(h,1);
plevels=NaN.*ones(h,1);
var_s_Thetalevels=[];
Thetalevels=[];


% exclude mixed layer that has not been mapped ---

jj=find(isnan(la_ptmp)==1);
PRES(jj)=NaN;
SAL(jj)=NaN;
PTMP(jj)=NaN;


% use only manually specified THETA & PRES range ---

if( isempty(use_theta_lt)==0 & isempty(use_theta_gt)==1 )
  jj=find(PTMP>use_theta_lt);
  PRES(jj)=NaN;
  SAL(jj)=NaN;
  PTMP(jj)=NaN;
end

if( isempty(use_theta_gt)==0 & isempty(use_theta_lt)==1 )
  jj=find(PTMP<use_theta_gt);
  PRES(jj)=NaN;
  SAL(jj)=NaN;
  PTMP(jj)=NaN;
end

if( isempty(use_theta_gt)==0 & isempty(use_theta_lt)==0 )
  if(use_theta_gt>use_theta_lt) %the middle band is excluded
    jj=find(PTMP<use_theta_gt&PTMP>use_theta_lt);
  else
    jj=find(PTMP<use_theta_gt|PTMP>use_theta_lt);
  end
  PRES(jj)=NaN;
  SAL(jj)=NaN;
  PTMP(jj)=NaN;
end

if( isempty(use_pres_lt)==0 & isempty(use_pres_gt)==1 )
  jj=find(PRES>use_pres_lt);
  PRES(jj)=NaN;
  SAL(jj)=NaN;
  PTMP(jj)=NaN;
end

if( isempty(use_pres_gt)==0 & isempty(use_pres_lt)==1 )
  jj=find(PRES<use_pres_gt);
  PRES(jj)=NaN;
  SAL(jj)=NaN;
  PTMP(jj)=NaN;
end

if( isempty(use_pres_gt)==0 & isempty(use_pres_lt)==0 )
  if(use_pres_gt>use_pres_lt) %he middle band is excluded
    jj=find(PRES<use_pres_gt&PRES>use_pres_lt);
  else
    jj=find(PRES<use_pres_gt|PRES>use_pres_lt);
  end
  PRES(jj)=NaN;
  SAL(jj)=NaN;
  PTMP(jj)=NaN;
end


% Define Thetalevels about 50 dbar apart ---
% Equally divide PRES into 50 dbar bins, then linearly interpolate PTMP onto
% these Preslevels. The Thetalevels are the mean of the interpolated PTMP.


%minTheta = ceil(min(min(PTMP))*10)/10;
%maxTheta = floor(max(max(PTMP))*10)/10;
minTheta = ceil(min(min(PTMP))*100)/100;  % GN
maxTheta = floor(max(max(PTMP))*100)/100;

if(minTheta<maxTheta) % no levels when no valid theta range --

increment=50;
maxpres=max(max(PRES));
minpres=min(min(PRES));
Preslevels=[minpres:increment:maxpres]';
if(length(Preslevels)<h) %if the usable water column is less than 500dbar thick,
  increment=floor((maxpres-minpres)/h); %hence fewer than 10 levels of 50dbar apart
  Preslevels=[minpres:increment:maxpres]';
end

interp_t=NaN.*ones(length(Preslevels),n);
for g=1:n
  ii=find(isnan(PRES(:,g))==0&isnan(PTMP(:,g))==0);
  for k=1:length(Preslevels)
     if( Preslevels(k)<max(PRES(ii,g)) & Preslevels(k)>min(PRES(ii,g)) )
       interp_t(k,g) = interp1( PRES(ii,g), PTMP(ii,g), Preslevels(k), 'linear');
     end
  end
end

Thetalevels = NaN*ones(length(Preslevels),1);
Thetalevel_indexes = NaN*ones(length(Preslevels),n);
for k=1:length(Preslevels)
  jj=find(isnan(interp_t(k,:))==0);
  if(isempty(jj)==0)
    Thetalevels(k)=mean(interp_t(k,jj));
  end
end


% find profile levels closest to theta levels ---
% At areas with temp inversions (e.g. Gulf of Alaska, Southern Ocean),
% PTMP is not unique, so this will pick out indexes from different depths for the same T,
% thus giving artificially high S var on T. This is ok since areas with
% temp inversions are also areas with naturally high S var on T, and therefore
% will not be selected as a min S var on T level anyway.

for i=1:n
  for j=1:length(Thetalevels)
   if(Thetalevels(j)<max(PTMP(:,i))&Thetalevels(j)>min(PTMP(:,i)))
     diffTheta = abs(PTMP(:,i)-Thetalevels(j));
     if isempty(find(~isnan(diffTheta)))
       Thetalevel_indexes(j,i) = NaN;
     else
       Thetalevel_indexes(j,i) = min(find(diffTheta==min(diffTheta)));
     end
   end
  end
end


% find S variance on theta levels ---

S_temp=NaN.*ones(length(Thetalevels),n);

for i=1:length(Thetalevels) % build the S matrix to find var & the max number of good profiles
  for j=1:n
   ti=Thetalevel_indexes(i,j);
   if ~isnan(ti)
     interval=max(ti-1,1):min(ti+1,m); %interval is one above and one below ti
     a=PTMP(ti,j) - PTMP(interval, j);
     if(PTMP(ti,j)>Thetalevels(i))
        gg=find(a>0);
        if( ~isempty(gg) )
           b=find(a==min(a(gg))); %find the level with min +ve diff
           ki=interval(b);
        else
           ki=ti;
        end
     end
     if(PTMP(ti,j)<Thetalevels(i))
        gg=find(a<0);
        if( ~isempty(gg) )
           b=find(-a==min(-a(gg))); %find the level with min -ve diff
           ki=interval(b);
        else
           ki=ti;
        end
     end
     if(PTMP(ti,j)==Thetalevels(i)) %this can happen when there is only 1 profile
        ki=ti;
     end
     if( ki~=ti&~isnan(SAL(ti,j))&~isnan(SAL(ki,j))&~isnan(PTMP(ti,j))&~isnan(PTMP(ki,j)) )
       S_temp(i,j) = interp1( [PTMP(ti,j), PTMP(ki,j)], [SAL(ti,j), SAL(ki,j)], Thetalevels(i) );
     else
       S_temp(i,j)=SAL(ti,j); % interpolate if possible because that is more accurate than using closest points
     end
   end
  end
end


numgood=NaN.*ones(length(Thetalevels),1);
percentSprofs=NaN.*ones(length(Thetalevels),1);
var_s_Thetalvls=NaN.*ones(length(Thetalevels),1);

for i=1:length(Thetalevels)
  good=find(isnan(S_temp(i,:))==0);
  numgood(i)=length(good);
  if(isempty(good)==0)
    var_s_Thetalvls(i) =  var(S_temp(i,good)); % only use S on theta level that has valid values
  end
end

for j=1:length(Thetalevels)
  if(max(numgood)~=0)
    percentSprofs(j) = numgood(j)/max(numgood); % maximum number of good data points on any theta level
  end
end
bad = find(percentSprofs<use_percent_gt);
var_s_Thetalvls(bad)=NaN;

var_s_Thetalevels=var_s_Thetalvls;


% select 10 thetas, these indexes will differ between profiles ---

for i=1:h
  if(isempty(find(isnan(var_s_Thetalvls)==0))==0)
    kk=min(find(var_s_Thetalvls==min(var_s_Thetalvls(find(isnan(var_s_Thetalvls)==0)) )));
    index(i,:) = Thetalevel_indexes(kk,:);
    tlevels(i) = Thetalevels(kk); %theta that has min S variance
    plevels(i) = Preslevels(kk); %corresponding nominal pressure
    var_s_Thetalvls(kk)=NaN; %block out this level --
  end
end


end % no levels when no valid theta range --


