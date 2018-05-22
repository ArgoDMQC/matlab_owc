
function plot_diagnostics_ow( pn_float_dir, pn_float_name, po_system_configuration )

%
% Annie Wong, 14 June 2011
% Breck Owens, October 2006
%--------------------------------------------------------------------------

%pn_float_dir='uw/';
%pn_float_name='R5902134';
%po_system_configuration = load_configuration( 'ow_config.txt' );


close all

% modify pn_float_name for title if the name contains '_' ------------

ii=find(pn_float_name=='_');
if(isempty(ii)==0)
  title_floatname = strcat( pn_float_name(1:ii-1), '\_', pn_float_name(ii+1:length(pn_float_name)));
else
  title_floatname = pn_float_name;
end


% load data from /float_source, /float_mapped, /float_calib, and config --------------

lo_float_source_data = load(fullfile( po_system_configuration.FLOAT_SOURCE_DIRECTORY, pn_float_dir, ...
  strcat( pn_float_name, po_system_configuration.FLOAT_SOURCE_POSTFIX ) ) );

PROFILE_NO = lo_float_source_data.PROFILE_NO;
LAT  = lo_float_source_data.LAT;
LONG = lo_float_source_data.LONG;
PRES = lo_float_source_data.PRES;
TEMP = lo_float_source_data.TEMP;
PTMP = lo_float_source_data.PTMP;
SAL  = lo_float_source_data.SAL;

if(isempty(find(isnan(PRES)==0))==0) % if no data exists, terminate here, no plots will be produced

lo_float_mapped_data = load( fullfile( po_system_configuration.FLOAT_MAPPED_DIRECTORY, pn_float_dir, ...
  strcat( po_system_configuration.FLOAT_MAPPED_PREFIX, pn_float_name, po_system_configuration.FLOAT_MAPPED_POSTFIX ) ) ) ;

mapped_sal  = lo_float_mapped_data.la_mapped_sal;
mapsalerrors = lo_float_mapped_data.la_mapsalerrors;
la_ptmp = lo_float_mapped_data.la_ptmp;
selected_hist =lo_float_mapped_data.selected_hist;

if(isempty(find(isnan(mapped_sal)==0))==0) % if mapping exists, terminate here, no plots will be produced

lo_float_calib_data = load( fullfile( po_system_configuration.FLOAT_CALIB_DIRECTORY, pn_float_dir, ...
  strcat( po_system_configuration.FLOAT_CALIB_PREFIX, pn_float_name, po_system_configuration.FLOAT_CALIB_POSTFIX ) ) );

cal_SAL = lo_float_calib_data.cal_SAL;
cal_SAL_err = lo_float_calib_data.cal_SAL_err;
pcond_factor = lo_float_calib_data.pcond_factor;
pcond_factor_err = lo_float_calib_data.pcond_factor_err;

lo_float_calseries = load( fullfile( po_system_configuration.FLOAT_CALIB_DIRECTORY, pn_float_dir, ...
  strcat( po_system_configuration.FLOAT_CALSERIES_PREFIX , pn_float_name, po_system_configuration.FLOAT_MAPPED_POSTFIX ) ) );

use_theta_gt = lo_float_calseries.use_theta_gt;
use_theta_lt = lo_float_calseries.use_theta_lt;
use_pres_gt = lo_float_calseries.use_pres_gt;
use_pres_lt = lo_float_calseries.use_pres_lt;
use_percent_gt = lo_float_calseries.use_percent_gt;

% load the station by station fits
load(fullfile( po_system_configuration.FLOAT_CALIB_DIRECTORY, pn_float_dir,...
   strcat( po_system_configuration.FLOAT_CALIB_PREFIX, pn_float_name, po_system_configuration.FLOAT_CALIB_POSTFIX ) ),'-regexp','^sta')


% plot the float locations (figure 1) -----------------------

load( fullfile( po_system_configuration.CONFIG_DIRECTORY, po_system_configuration.CONFIG_COASTLINES ), 'coastdata_x', 'coastdata_y' );

[m,n] = size(PRES);

figure
set(gcf,'defaultaxeslinewidth',2)
set(gcf,'defaultlinelinewidth',2)
set(gcf,'defaultaxesfontsize',16)

colormap(jet(n));
c=colormap;

x=[];
y=[];
if(isempty(selected_hist)==0)
  x=selected_hist(:,1);
  y=selected_hist(:,2);
end

if( max(LONG)+30>360 | min(LONG)-30<0 ) % if there are float data within 30 deg lon of the prime meridian
  if( isempty( find(LAT<-50&LONG>250&LONG<285) )==0 )
     jj=find(LONG<LONG(1)-30); % and if there are float data west of the Drake Passage
     LONG(jj)=LONG(jj)+360;    % then use the first float location point minus 30 deg lon
     kk=find(x<LONG(1)-30);    % as the western edge of the map
     x(kk)=x(kk)+360;
     mm=find(coastdata_x<LONG(1)-30);
     coastdata_x(mm)=coastdata_x(mm)+360;
     nn=find(coastdata_x<LONG(1)-30+1|coastdata_x>LONG(1)-30+360-1); %remove coastline wraparound
     coastdata_x(nn)=NaN;
  else
     jj=find(LONG<250);     % or else if there are no float data west of the Drake Passage
     LONG(jj)=LONG(jj)+360; % then use 250 deg lon as the western edge of the map
     kk=find(x<250);
     x(kk)=x(kk)+360;
     mm=find(coastdata_x<250);
     coastdata_x(mm)=coastdata_x(mm)+360;
     nn=find(coastdata_x<251|coastdata_x>609); %remove coastline wraparound
     coastdata_x(nn)=NaN;
  end
end


plot(LONG,LAT,'r-');
hold on
plot(x,y,'b.')
legend('float','historical points','Location','Best')
plot(coastdata_x,coastdata_y,'k.-');

for i=1:n
  h=plot(LONG(i),LAT(i),'+');
  set(h,'color',c(i,:));
  j=text(LONG(i),LAT(i),int2str(PROFILE_NO(i)));
  set(j,'color',c(i,:),'fontsize',12,'hor','cen');
end
axis([min(LONG)-30, max(LONG)+30, min(LAT)-20, max(LAT)+20])
set(gca,'FontSize',12)
xlabel('Longitude');
ylabel('Latitude');
title( strcat(title_floatname, ' profile locations with historical data' ) );

h=gca;
xticks=get(h,'XTick');
xticklabels=get(h,'XTickLabel');
ii=find(xticks>360);
xticks(ii)=xticks(ii)-360;
for j=1:length(ii)
  enough=num2str(ones(3,1));
  c=char(num2str(xticks(ii(j))),enough');
  xticklabels(ii(j),:)=c(1,:);
end
set(gca,'XTickLabel',xticklabels);

drawnow
set(gcf,'papertype','usletter','paperunits','inches','paperorientation','portrait','paperposition',[.25,.75,8,9.5]);
print('-depsc ', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, pn_float_dir, pn_float_name, '_1.eps'));


% plot the uncalibrated theta-S curves from the float (figure 2) --------

figure;
set(gcf,'defaultaxeslinewidth',2)
set(gcf,'defaultlinelinewidth',2)
set(gcf,'defaultaxesfontsize',14)

jj=[1:ceil(n/30):n]; % the legend can only fit 30 profiles
ok = [];
for k=1:length(jj)
  if ~isempty(find(isfinite(SAL(:,jj(k)))))
      ok = [ok jj(k)];
  end
end
jj = ok;

colormap(jet(n));
c=colormap;
qq = plot(PTMP(:,jj), SAL(:,jj));
for i=1:length(jj)
  set(qq(i),'color',c(jj(i),:)) ;
end
hold on
legend(int2str([PROFILE_NO(jj)]'), 'Location', 'NorthEastOutside')

for i=1:n % plot all remaining profiles
  qq = plot( PTMP(:,i), SAL(:,i) );
  set(qq,'color',c(i,:)) ;
end

[tlevels, plevels, index, var_s_Thetalevels, Thetalevels] = find_10thetas( SAL, PTMP, PRES, la_ptmp, use_theta_gt, use_theta_lt, use_pres_gt, use_pres_lt, use_percent_gt);

for i=1:n
  b = find( isnan(index(:,i))==0 );
  a = index(b,i);
  if ~isempty(a)
   h = errorbar( la_ptmp(a,i), mapped_sal(a,i), mapsalerrors(a,i), 'o', 'color', c(i,:) );
  end
end

view([90 -90]);
set(gca,'FontSize',12)
ylabel('Salinity (PSS-78)');
xlabel('\theta ^{\circ} C');

max_s=max([max(SAL),max(mapped_sal)])+.1;
min_s=min([min(SAL),min(mapped_sal)])-.1;
max_t=max(max(PTMP))+1;
min_t=min(min(PTMP))-1;
axis([min_t,max_t,min_s,max_s]);

drawnow
title( strcat( title_floatname, ' uncalibrated float data (-) and mapped salinity (o) with objective errors' ) );
set(gcf,'papertype','usletter','paperunits','inches','paperorientation','portrait','paperposition',[.25,.75,8,9.5]);
print('-depsc ', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, pn_float_dir, pn_float_name, '_2.eps'));


% calibration curve (figure 3) --------------------------

if(isempty(find(isnan(cal_SAL)==0))==0) % if no cal exists, terminate here, no plot will be produced

Soffset=cal_SAL-SAL;
avg_Soffset=NaN.*ones(1,n);
avg_Soffset_err=NaN.*ones(1,n);
Staoffset=sta_SAL-SAL;
avg_Staoffset=NaN.*ones(1,n);
avg_Staoffset_err=NaN.*ones(1,n);

for i=1:n
   ii=[];
   ii=find(isnan(Soffset(:,i))==0);
   if ~isempty(ii)
       avg_Soffset(i)=mean(Soffset(ii,i));
       avg_Soffset_err(i)=mean(cal_SAL_err(ii,i));
   else
       avg_Soffset(i) = NaN;
       avg_Soffset_err(i) = NaN;
   end
   ii=find(isnan(Staoffset(:,i))==0);
   if ~isempty(ii)
       avg_Staoffset(i)=mean(Staoffset(ii,i));
       avg_Staoffset_err(i)=mean(sta_SAL_err(ii,i));
   else
       avg_Staoffset(i) = NaN;
       avg_Staoffset_err(i) = NaN;
   end
end

figure
set(gcf,'defaultaxeslinewidth',2)
set(gcf,'defaultlinelinewidth',2)
set(gcf,'defaultaxesfontsize',16)
subplot(2,1,1)
plot(PROFILE_NO, pcond_factor, 'b-');
hold on
plot(PROFILE_NO, pcond_factor, 'g-');
% plot station by station fit
ok = find(isfinite(sta_mean));
plot(PROFILE_NO(ok), sta_mean(ok), 'r-');
legend('2 x cal error','1 x cal error','1-1 profile fit', 'Location', 'Best');
errorbar(PROFILE_NO, pcond_factor, 2*pcond_factor_err,'b')
errorbar(PROFILE_NO, pcond_factor, pcond_factor_err,'g*-')
plot(PROFILE_NO(ok), sta_mean(ok), 'r-');

plot( [0, max(PROFILE_NO)+1], [1,1], 'k-')
axis([ 0, max(PROFILE_NO)+1, min([pcond_factor-pcond_factor_err,1])-.0004, max([pcond_factor+pcond_factor_err,1])+.0004 ])
set(gca,'FontSize',12)
ylabel('r') % multiplicative term has no units
title( strcat(title_floatname, ' potential conductivity (mmho/cm) multiplicative correction r with errors') );

subplot(2,1,2)
plot(PROFILE_NO, avg_Soffset, 'b-');
hold on
plot(PROFILE_NO, avg_Soffset, 'g-');
% Plot station by station fit
ok = find(isfinite(avg_Staoffset));
plot(PROFILE_NO(ok), avg_Staoffset(ok), 'r-');
legend('2 x cal error','1 x cal error','1-1 profile fit', 'Location', 'Best');
errorbar(PROFILE_NO, avg_Soffset, 2*avg_Soffset_err,'b')
errorbar(PROFILE_NO, avg_Soffset, avg_Soffset_err,'go-')
plot(PROFILE_NO(ok), avg_Staoffset(ok), 'r-');

axis([ 0, max(PROFILE_NO)+1, min([avg_Soffset-avg_Soffset_err,0])-.02, max([avg_Soffset+avg_Soffset_err,0])+.02 ])
plot( [0, max(PROFILE_NO)+1], [0,0], 'k-')
set(gca,'FontSize',12)
xlabel('float profile number');
ylabel('\Delta S')
title( strcat(title_floatname, ' vertically-averaged salinity (PSS-78) additive correction \Delta S with errors') );

drawnow
set(gcf,'papertype','usletter','paperunits','inches','paperorientation','portrait','paperposition',[.25,.75,8,9.5]);
print('-depsc ', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, pn_float_dir, pn_float_name, '_3.eps'));


% plot the calibrated theta-S curves from the float (figure 4) --------------------------

figure;
set(gcf,'defaultaxeslinewidth',2)
set(gcf,'defaultlinelinewidth',2)
set(gcf,'defaultaxesfontsize',14)

jj=[1:ceil(n/30):n]; % the legend can only fit 30 profiles
ok = [];% check to make sure we have choosen profiles with good data
for k=1:length(jj)
  if ~isempty(find(isfinite(SAL(:,jj(k)))))
      ok = [ok jj(k)];
  end
end
jj = ok;

colormap(jet(n));
c=colormap;
qq=plot( PTMP(:,jj), cal_SAL(:,jj) );
for i=1:length(jj)
 set(qq(i),'color',c(jj(i),:)) ;
end
hold on;
legend(int2str([PROFILE_NO(jj)]'), 'Location', 'NorthEastOutside')

for i=1:n % plot all remaining profiles
  qq = plot( PTMP(:,i), cal_SAL(:,i) );
  set(qq,'color',c(i,:)) ;
end

for i=1:n
  b = find( isnan(index(:,i))==0 );
  a = index(b,i);
  if ~isempty(a)
    h = errorbar( la_ptmp(a,i), mapped_sal(a,i), mapsalerrors(a,i), 'o', 'color', c(i,:));
  end
end

view([90 -90]);
set(gca,'FontSize',12)
xlabel('\theta ^{\circ} C')
ylabel('Salinity (PSS-78)')

max_s=max([max(max(cal_SAL)),max(max(mapped_sal))])+.1;
min_s=min([min(min(cal_SAL)),min(min(mapped_sal))])-.1;
max_t=max(max(PTMP))+1;
min_t=min(min(PTMP))-1;
if(isnan(min_s)==0)
  axis([min_t,max_t,min_s,max_s]);
end

drawnow
title( strcat(title_floatname, ' calibrated float data (-) and mapped salinity (o) with objective errors' ) );
set(gcf,'papertype','usletter','paperunits','inches','paperorientation','portrait','paperposition',[.25,.75,8,9.5]);
print('-depsc ', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, pn_float_dir, pn_float_name, '_4.eps'));


% Brian King's plot: salinity anomaly time series on theta levels (figure 5) ------------

figure
set(gcf,'defaultaxeslinewidth',2)
set(gcf,'defaultlinelinewidth',2)
set(gcf,'defaultaxesfontsize',16)

fl.useqc = '0';
fl.plot = 1;
fl.yaxes = [2 5 20];
d.PSAL = SAL;
d.TEMP = TEMP;
d.PRES = PRES;
d.PSAL_QC = zeros(m,n);
d.TEMP_QC = zeros(m,n);
d.PRES_QC = zeros(m,n);
d.LONGITUDE = LONG;
d.LATITUDE = LAT;
d.PROFILE_NO = PROFILE_NO;
fl = anom(d,fl); % Brian King's routine
subplot('position',[.1 .45 .8 .35])
title(['       Salinity anom on theta.    ' title_floatname])

drawnow
set(gcf,'papertype','usletter','paperunits','inches','paperorientation','portrait','paperposition',[.25,.5,8,10]);
print('-depsc ', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, pn_float_dir, pn_float_name, '_5.eps'));


% plot salinity time series on theta levels with the smallest S variance (figure 6) ------------

[tlevels, plevels, index, var_s_Thetalevels, Thetalevels] = find_10thetas( SAL, PTMP, PRES, la_ptmp, use_theta_gt, use_theta_lt, use_pres_gt, use_pres_lt, use_percent_gt);

tplot=[1:2];
p=length(tplot);
Sint=NaN.*ones(p,n);
Smap=NaN.*ones(p,n);
Smaperr=NaN.*ones(p,n);
Scal=NaN.*ones(p,n);
Scalerr=NaN.*ones(p,n);
Thetalevel_indexes=NaN.*ones(p,n);

trimPRES=PRES;  % use only manually specified THETA & PRES range ---
trimSAL=SAL;
trimPTMP=PTMP;
trim_mapped_sal=mapped_sal;
trim_mapsalerrors=mapsalerrors;
trim_cal_SAL=cal_SAL;
trim_cal_SAL_err=cal_SAL_err;

jj=find(isnan(la_ptmp)==1);
trimPRES(jj)=NaN;
trimSAL(jj)=NaN;
trimPTMP(jj)=NaN;
trim_mapped_sal(jj)=NaN;
trim_mapsalerrors(jj)=NaN;
trim_cal_SAL(jj)=NaN;
trim_cal_SAL_err(jj)=NaN;

if( isempty(use_theta_lt)==0 & isempty(use_theta_gt)==1 )
  jj=find(trimPTMP>use_theta_lt);
  trimPRES(jj)=NaN;
  trimSAL(jj)=NaN;
  trimPTMP(jj)=NaN;
  trim_mapped_sal(jj)=NaN;
  trim_mapsalerrors(jj)=NaN;
  trim_cal_SAL(jj)=NaN;
  trim_cal_SAL_err(jj)=NaN;
end

if( isempty(use_theta_gt)==0 & isempty(use_theta_lt)==1 )
  jj=find(trimPTMP<use_theta_gt);
  trimPRES(jj)=NaN;
  trimSAL(jj)=NaN;
  trimPTMP(jj)=NaN;
  trim_mapped_sal(jj)=NaN;
  trim_mapsalerrors(jj)=NaN;
  trim_cal_SAL(jj)=NaN;
  trim_cal_SAL_err(jj)=NaN;
end

if( isempty(use_theta_gt)==0 & isempty(use_theta_lt)==0 )
  if(use_theta_gt>use_theta_lt) %the middle band is excluded
    jj=find(trimPTMP<use_theta_gt&trimPTMP>use_theta_lt);
  else
    jj=find(trimPTMP<use_theta_gt|trimPTMP>use_theta_lt);
  end
  trimPRES(jj)=NaN;
  trimSAL(jj)=NaN;
  trimPTMP(jj)=NaN;
  trim_mapped_sal(jj)=NaN;
  trim_mapsalerrors(jj)=NaN;
  trim_cal_SAL(jj)=NaN;
  trim_cal_SAL_err(jj)=NaN;
end

if( isempty(use_pres_lt)==0 & isempty(use_pres_gt)==1 )
  jj=find(trimPRES>use_pres_lt);
  trimPRES(jj)=NaN;
  trimSAL(jj)=NaN;
  trimPTMP(jj)=NaN;
  trim_mapped_sal(jj)=NaN;
  trim_mapsalerrors(jj)=NaN;
  trim_cal_SAL(jj)=NaN;
  trim_cal_SAL_err(jj)=NaN;
end

if( isempty(use_pres_gt)==0 & isempty(use_pres_lt)==1 )
  jj=find(trimPRES<use_pres_gt);
  trimPRES(jj)=NaN;
  trimSAL(jj)=NaN;
  trimPTMP(jj)=NaN;
  trim_mapped_sal(jj)=NaN;
  trim_mapsalerrors(jj)=NaN;
  trim_cal_SAL(jj)=NaN;
  trim_cal_SAL_err(jj)=NaN;
end

if( isempty(use_pres_gt)==0 & isempty(use_pres_lt)==0 )
  if(use_pres_gt>use_pres_lt) %he middle band is excluded
    jj=find(trimPRES<use_pres_gt&trimPRES>use_pres_lt);
  else
    jj=find(trimPRES<use_pres_gt|trimPRES>use_pres_lt);
  end
  trimPRES(jj)=NaN;
  trimSAL(jj)=NaN;
  trimPTMP(jj)=NaN;
  trim_mapped_sal(jj)=NaN;
  trim_mapsalerrors(jj)=NaN;
  trim_cal_SAL(jj)=NaN;
  trim_cal_SAL_err(jj)=NaN;
end

for i=1:n
  for j=tplot
   if(tlevels(j)<max(trimPTMP(:,i))&tlevels(j)>min(trimPTMP(:,i)))
    diffTheta = abs(trimPTMP(:,i)-tlevels(j));
    if isempty(find(~isnan(diffTheta)))
      Thetalevel_indexes(j,i) = NaN;
    else
      Thetalevel_indexes(j,i) = min(find(diffTheta==min(diffTheta)));
    end
   end
  end
end

for i=tplot % build the S matrix for plotting
  for j=1:n
   ti=Thetalevel_indexes(i,j);
   if ~isnan(ti)
     interval=max(ti-1,1):min(ti+1,m); %interval is one above and one below ti
     a = trimPTMP(ti,j) - trimPTMP(interval, j);
     if( trimPTMP(ti,j)>tlevels(i) )
        gg=find(a>0);
        if( ~isempty(gg) )
           b=find(a==min(a(gg))); %find the level with min +ve diff
           ki=interval(b);
        else
           ki=ti;
        end
     end
     if( trimPTMP(ti,j)<tlevels(i) )
        gg=find(a<0);
        if( ~isempty(gg) )
           b=find(-a==min(-a(gg))); %find the level with min -ve diff
           ki=interval(b);
        else
           ki=ti;
        end
     end
     if( trimPTMP(ti,j)==tlevels(i) )
        ki=ti;
     end
     if( ki~=ti&~isnan(trimSAL(ti,j))&~isnan(trimSAL(ki,j))&~isnan(trimPTMP(ti,j))&~isnan(trimPTMP(ki,j)) )
       Sint(i,j) = interp1( [trimPTMP(ti,j), trimPTMP(ki,j)], [trimSAL(ti,j), trimSAL(ki,j)], tlevels(i) );
     else
       Sint(i,j) = trimSAL(ti,j); % interpolate if possible because that is more accurate than using closest points
     end
     if( ki~=ti&~isnan(trim_mapped_sal(ti,j))&~isnan(trim_mapped_sal(ki,j))&~isnan(trimPTMP(ti,j))&~isnan(trimPTMP(ki,j)) )
       Smap(i,j) = interp1( [trimPTMP(ti,j), trimPTMP(ki,j)], [trim_mapped_sal(ti,j), trim_mapped_sal(ki,j)], tlevels(i) );
       Smaperr(i,j) = interp1( [trimPTMP(ti,j), trimPTMP(ki,j)], [trim_mapsalerrors(ti,j), trim_mapsalerrors(ki,j)], tlevels(i) );
     else
       Smap(i,j)=trim_mapped_sal(ti,j); % interpolate if possible because that is more accurate than using closest points
       Smaperr(i,j)=trim_mapsalerrors(ti,j); % interpolate if possible because that is more accurate than using closest points
     end
     if( ki~=ti&~isnan(trim_cal_SAL(ti,j))&~isnan(trim_cal_SAL(ki,j))&~isnan(trimPTMP(ti,j))&~isnan(trimPTMP(ki,j)) )
       Scal(i,j) = interp1( [trimPTMP(ti,j), trimPTMP(ki,j)], [trim_cal_SAL(ti,j), trim_cal_SAL(ki,j)], tlevels(i) );
       Scalerr(i,j) = interp1( [trimPTMP(ti,j), trimPTMP(ki,j)], [trim_cal_SAL_err(ti,j), trim_cal_SAL_err(ki,j)], tlevels(i) );
     else
       Scal(i,j)=trim_cal_SAL(ti,j); % interpolate if possible because that is more accurate than using closest points
       Scalerr(i,j)=trim_cal_SAL_err(ti,j); % interpolate if possible because that is more accurate than using closest points
     end
   end
  end
end

figure
set(gcf,'defaultaxeslinewidth',2)
set(gcf,'defaultlinelinewidth',2)
set(gcf,'defaultaxesfontsize',16)

for k=1:p
j=tplot(k);
subplot(p,1,k)
if(isempty(Sint)==0)
  plot(PROFILE_NO,Sint(j,:),'b*-');
  hold on
  plot(PROFILE_NO,Smap(j,:),'r');
  plot(PROFILE_NO,Scal(j,:),'g');
  mm=find(isfinite(Scal(j,:))==1); ll=PROFILE_NO; ll=ll(mm); kk=Scal(j,mm); nn=Scalerr(j,mm);
  h=fill([ll,fliplr(ll)],[kk+nn,fliplr([kk-nn])],'g');
  set(h,'EdgeColor','g');
  errorbar(PROFILE_NO,Smap(j,:),Smaperr(j,:),'r-')
  plot(PROFILE_NO,Sint(j,:),'b*-');
  SMIN = min([Sint(j,:),Scal(j,:),Smap(j,:)]);
  SMAX = max([Sint(j,:),Scal(j,:),Smap(j,:)]);
  if isfinite(SMIN) & isfinite(SMAX)
    axis([0, max(PROFILE_NO)+1, SMIN-.05, SMAX+.05 ])
  end
  set(gca,'FontSize',12)
  title( strcat(title_floatname, ' salinities with error on \theta= ', num2str(tlevels(j)), '^{\circ}C' ) );
  ylabel('PSS-78')
end
end
legend('uncal float','mapped salinity','cal float w/1xerr.', 'Location', 'Best');
xlabel('float profile number');

drawnow
set(gcf,'papertype','usletter','paperunits','inches','paperorientation','portrait','paperposition',[.25,.75,8,9.5]);
print('-depsc ', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, pn_float_dir, pn_float_name, '_6.eps'));


% Brian King's plot: salinity anomaly time series on theta levels (figure 7) ------------

figure
set(gcf,'defaultaxeslinewidth',2)
set(gcf,'defaultlinelinewidth',2)
set(gcf,'defaultaxesfontsize',16)

fl.useqc = '0';
fl.plot = 1;
fl.yaxes = [2 5 20];
d.PSAL = cal_SAL;
d.TEMP = TEMP;
d.PRES = PRES;
d.PSAL_QC = zeros(m,n);
d.TEMP_QC = zeros(m,n);
d.PRES_QC = zeros(m,n);
d.LONGITUDE = LONG;
d.LATITUDE = LAT;
d.PROFILE_NO = PROFILE_NO;
fl = anom(d,fl); % Brian King's routine
subplot('position',[.1 .45 .8 .35])
title(['Calibrated salinity anom on theta. ' title_floatname])

drawnow
set(gcf,'papertype','usletter','paperunits','inches','paperorientation','portrait','paperposition',[.25,.5,8,10]);
print('-depsc ', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, pn_float_dir, pn_float_name, '_7.eps'));


% Paul Robbins' analyse variance plot (figure 8) ------------

figure
set(gcf,'defaultaxeslinewidth',1)
set(gcf,'defaultlinelinewidth',1)
set(gcf,'defaultaxesfontsize',12)

% plot t-s profile
subplot(222)
plot(SAL,PTMP,'b-');
x = get(gca,'xlim');
y = get(gca,'ylim');
xlabel('PSS-78')
title(strcat('OW chosen levels - ', pn_float_name));
for i=1:10
  hold on
  plot(x,[tlevels(i) tlevels(i)] ,'g-');
end

% plot s var on t
subplot(221)
plot(var_s_Thetalevels,Thetalevels,'b-')
x = get(gca,'xlim');
for i=1:10
  hold on
  plot(x,[tlevels(i) tlevels(i)] ,'g-');
end
title('Salinity Variance on Theta')
ylabel('Potential temp (^{\circ}C)')
xlabel('salinity variance')
%set(gca,'ylim',y)

% plot p-t profile
subplot(223)
plot(PTMP,-PRES,'b-');
x = get(gca,'xlim');
xlabel('^{\circ}C')
ylabel('Pressure (dbar)')
title(strcat('OW chosen levels - ', pn_float_name));
for i=1:10
  hold on
  plot(x,[-plevels(i) -plevels(i)] ,'g-');
end

% plot p-s profile
subplot(224)
plot(SAL,-PRES,'b-');
x = get(gca,'xlim');
xlabel('PSS-78')
title(strcat('OW chosen levels - ', pn_float_name));
for i=1:10
  hold on
  plot(x,[-plevels(i) -plevels(i)] ,'g-');
end


drawnow
set(gcf,'papertype','usletter','paperunits','inches','paperorientation','portrait','paperposition',[.5,.25,8,10.25]);
print('-depsc ', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, pn_float_dir, pn_float_name, '_8.eps'));


end %if(isempty(find(isnan(cal_SAL)==0))==0) ---------------

end %if(isempty(find(isnan(mapped_sal)==0))==0) -----------

end %if(isempty(find(isnan(PRES)==0))==0) ------------------


