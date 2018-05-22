
function fl = anom(d,fl)

% B. King, April 2005.
%
% Brian King's routine for plotting salinity anomaly time series.
% Needs 3 other routines: axmerc.m, setdefault.m, make_colorbar.m.
%
% assume all profiles for a float have been read in from a single netcdf file
% the argument is the structure array created when reading the netcdf file
%
% use:
% fl = anom(d)   or   
%
% fl.useqc = '0';
% fl.plot = 0;   to switch off plotting
% fl.plot = 1; to switch on plotting
% fl.yaxes = [2 5 20]; define the theta levels for plot limits and to divide panels
% fl = anom(d,fl)
%
% The structure d must have the following fields
% d.PLATFORM_NUMBER
% d.PSAL
% d.TEMP
% d.PRES
% d.PSAL_QC
% d.TEMP_QC
% d.PRES_QC
%
% The structure fl is used to pass back some of the working variables as follows:
%
%  fl.theta_base;  theta levels  300x1
%  fl.Smean;  S reference profile 300x1
%  fl.Sanom;  S anomaly array     300xn
%  fl.Sint;   S interpolated array 300xn
%
% the flag fl.useqc is a string variable. The contens of the string are
% the QC flags that are to be retained. eg fl.useqc = '1' retains only data
% with QC flag '1' and sets the rest to NaN. This is the default.
% Use fl.useqc = '0' to keep all data.
%
% note that the arrays contain a large number of zeros
% for depth levels greater than 100. These can be sorted out with
% something like snan = find(S < 1) and T(snan) = nan; P(snan) = nan; S(snan) = nan; etc


theta_base = 0.1:.1:30;
fl.theta_base = theta_base';
%AW plat = d.PLATFORM_NUMBER;
%AW platchar = plat(1,:);
%AW platchar(findstr(platchar,' '))=[];
clear S T P TH
S = d.PSAL;
T = d.TEMP;
P = d.PRES;
PROFILE_NO=d.PROFILE_NO; %AW
%S = S';
%T = T';
%P = P';
%S(find(S < 30)) = nan;

fl = setdefault(fl,'useqc','1'); %fl.useqc = '1' by default. ie good data only
fl = setdefault(fl,'plot',1);    %fl.plot = 1 by default. ie plotting on.
fl = setdefault(fl,'yaxes',[2 5 20]);

%AW disp(['platform = ' platchar '   useqc = ' num2str(fl.useqc) '   plot = ' num2str(fl.plot)])
%AWdisp(['y range for two panels of plot are: [' num2str(fl.yaxes(1:2)) '] [' num2str(fl.yaxes(2:3)) ']']);

%use qc flags if required; this is default unless fl.useqc passed into function
if strcmp(fl.useqc,'0')
    %retain all data; no qc flag processing
else
    sqc = d.PSAL_QC; sqc = sqc';
    tqc = d.TEMP_QC; tqc = tqc';
    pqc = d.PRES_QC; pqc = pqc';
    for row = 1:size(sqc,1)
        srownan = nan*ones(1,size(sqc,2));
        trownan = srownan;
        prownan = srownan;
        for flag = 1:length(fl.useqc) 
            %for each element of fl.useqc, find the matching flags to keep
            skeep = findstr(sqc(row,:),fl.useqc(flag));
            srownan(skeep) = 0;
            tkeep = findstr(tqc(row,:),fl.useqc(flag));
            trownan(tkeep) = 0;
            pkeep = findstr(pqc(row,:),fl.useqc(flag));
            prownan(pkeep) = 0;
        end 
        S(row,:) = S(row,:)+srownan; %discard unwanted data by adding 0 (keep) or nan (discard)
        T(row,:) = T(row,:)+trownan;
        P(row,:) = P(row,:)+prownan;
    end
end

TH = sw_ptmp(S,T,P,0);
Sint = nan*ones(length(theta_base),size(S,2));
Smean = squeeze(Sint(:,1));
Sstd = Smean;
Sanom = Sint;
prof_range = nan*ones(2,1);

    
    [i,j] = find(~isnan(TH(:,:)));
    prof_range(1) = min(j);
    prof_range(2) = max(j);    %This is the prof range that contains good data.
    
    for k = prof_range(1):prof_range(2)
        th = TH(:,k);
        s = S(:,k);
        w1 = [th s];
        w2 = sortrows(w1,1);
        [xx,xi,xj] = unique(w2(:,1)); % I once found a repeated level in a PROVOR profile
        clear w3
        w3(:,1) = w2(xi,1);
        w3(:,2) = w2(xi,2);
        w2 = w3;
        igood = find(~isnan(w2(:,1)) & ~isnan(w2(:,2)));
        if(length(igood) < 3)
            %insist on at least 3 good levels in the profile !
        else
            sint = interp1(w2(igood,1),w2(igood,2),theta_base(:));
            Sint(:,k) = sint;
        end
    end
%AW    disp('End of interpolation')
    
    s1 = Sint';
%AW    smed = nanmedian(s1);
%AW    sstd = nanstd(s1);

%--added by AW for people who don't have nanmedian & nanstd ---
    [m,n]=size(s1);
    smed=NaN.*ones(1,n);
    sstd=NaN.*ones(1,n);
    for i=1:n
      jj=find(isnan(s1(:,i))==0);
      if(isempty(jj)==0)
        smed(i)=median( s1(jj,i) );
        sstd(i)=std( s1(jj,i) );
      end
    end
%-----------------------------------------------------AW

    Sstd = sstd';
    Smean = smed';
    
   for k = prof_range(1):prof_range(2)
       Sanom(:,k) = Sint(:,k)-smed';
   end
   fl.Sanom = Sanom;
   fl.Sint = Sint;
   fl.Smean = Smean;

%AW disp('End of calculating anomaly')

if(fl.plot ~= 1) 
    return 
end
%AW disp('Plotting first contour section')

%make_colorbar


colorbar_limits = [-.1 .1];
color_level_boundaries = [ -.06 -.04 -.02 -.01 -.005 .005 .01 .02 .04 .06];
colortable = [
         0         0    0.7000
         0         0    0.9375
         0    0.3125    1.0000
         0    0.6875    1.0000
    0.8000    1.0000    1.0000
    1.0000    1.0000    1.0000
    1.0000    1.0000    0.6000
    1.0000    0.8750         0
    1.0000    0.6000         0
    1.0000    0.4000         0
    1.0000    0.1250         0
    ];
%define one more color than there are coloe level boundaries
cmap = make_colorbar(colorbar_limits,color_level_boundaries,colortable);

clev = [-.1 -.08 -.06 -.04 -.02 -.01 -.005 0 .005 .01 .02 .04 .06 .08 .1]; %contour levels

%    figure
    set(gcf,'defaultaxeslinewidth',2);
   % set(gcf,'defaultlinelinewidth',2)
    set(gcf,'defaultaxesfontsize',16);
    orient tall
    x1 = prof_range(1);
    x2 = prof_range(2);
    xx = PROFILE_NO(x1:x2); %aw
    yrange = fl.yaxes; %yrange defines theta levels for bottom, middle and top of plot
    y1 = find(yrange(1) < theta_base & theta_base < yrange(2));
    y2 = find(yrange(2) < theta_base & theta_base < yrange(3));
    
    subplot('position',[.1 .45 .8 .35])
    splot = Sanom(y2,x1:x2);
    smin = min(min(splot));
    smax = max(max(splot));
    if ~isnan(smin)
%aw        [c,h] = contourf(x1:x2,theta_base(y2),Sanom(y2,x1:x2),clev,'k-');
%breck--
        if length(y2) > 1 & length([x1:x2]) > 1
            [c,h] = contourf(xx,theta_base(y2),Sanom(y2,x1:x2),clev,'k-');
        else
            plot(nan,nan);
        end
    else
        plot(nan,nan) % if no data to plot, then plot something, so that a set of axes are drawn
    end
    %bakmap2 = load('/local/users/pstar/data/other_software/bak_NCDF/bak_scripts/bakmap2');
    colormap(cmap);
    caxis(colorbar_limits)
    %colorbar('horiz')
    ax = axis;
%aw    xaxes(1) = 10*floor(x1/10);
%aw    xaxes(2) = 10*ceil(x2/10);
    xaxes(1) = 10*floor(xx(1)/10);
    xaxes(2) = 10*ceil(xx(length(xx))/10);
    axis([xaxes yrange(2:3)]);
%AW    title(['       Salinity anom on theta.    WMO ' platchar])
    hold on
    if ~isnan(smin)
%aw        [c,hline] = contour(xx1:xx2,theta_base(y2),Sanom(y2,x1:x2),[0  0],'k-');
%breck--
        if length(y2) > 1 & length([x1:x2]) > 1
            [c,hline] = contour(xx,theta_base(y2),Sanom(y2,x1:x2),[0  0],'k-');
            set(hline,'linewidth',3)
        end
    end
    %xlabel('profile number')
    %ylabel(['theta'])
    
%AW    disp('Plotting second contour section')
 
    subplot('position',[.1 .1 .8 .35])
    
    splot = Sanom(y1,x1:x2);
    smin = min(min(splot));
    smax = max(max(splot));
    if ~isnan(smin)
%aw        [c,h] = contourf(xx1:xx2,theta_base(y1),Sanom(y1,x1:x2),clev,'k-');
%breck--
        if length(y1) > 1 & length([x1:x2]) > 1
            [c,h] = contourf(xx,theta_base(y1),Sanom(y1,x1:x2),clev,'k-');
        end
    else
        plot(nan,nan) %plot something, so that a set of axes are drawn
    end
    %bakmap2 = load('/local/users/pstar/data/other_software/bak_NCDF/bak_scripts/bakmap2');
    colormap(cmap);
    caxis(colorbar_limits)
    colorbar('horiz')
    ax = axis;
    axis([xaxes yrange(1:2)]);
   % title(['Salinity anom on theta. Index: ' sprintf('%3d',index) '  WMO ' sprintf('%07d',float_list(index))])
    hold on
    if ~isnan(smin)
%aw        [c,hline] = contour(xx1:xx2,theta_base(y1),Sanom(y1,x1:x2),[0  0],'k-');
%breck--
        if length(y1) > 1 & length([x1:x2]) > 1
            [c,hline] = contour(xx,theta_base(y1),Sanom(y1,x1:x2),[0  0],'k-');
            set(hline,'linewidth',3)
        end
    end    
    xlabel('profile number')
    ylabel(['theta'])
   
%AW    disp('Plotting track')
    subplot('position',[.1 .83 .15 .15])
    lon = d.LONGITUDE;
    lat = d.LATITUDE;
    plot(lon(x1:x2),lat(x1:x2),'+-')
    hold on
    grid on
    axmerc
    plot(lon(x1),lat(x1),'r*')
    plot(lon(x2),lat(x2),'ko')
%AW    psname = ['sanom_' platchar '.ps'];
    %eval(['print -dpsc ' psname])
 %{   
    figure
    %track plot
    set(gcf,'defaultaxeslinewidth',2)
    set(gcf,'defaultlinelinewidth',2)
    set(gcf,'defaultaxesfontsize',16)
    disp('track plot')
    lon = d.LONGITUDE;
    lat = d.LATITUDE;
    plot(lon(x1:x2),lat(x1:x2),'+-')
    hold on
    grid on
    axmerc    
%AW    title(['  WMO ' platchar])
    plot(lon(x1),lat(x1),'r*')
    plot(lon(x2),lat(x2),'ko')
    xlabel('Longitude');
    ylabel('Latitude');
%AW    psname = ['pos_' platchar '.ps'];
    %eval(['print -dpsc ' psname])
    
    figure
    %theta-S plot
    set(gcf,'defaultaxeslinewidth',2)
    set(gcf,'defaultlinelinewidth',2)
    set(gcf,'defaultaxesfontsize',16)
    disp('theta-S plot')
    plot(Sint,theta_base)
    hold on
    grid on
  %  fl.Sint = Sint;
  %  fl.theta_base = theta_base;
%AW    title(['  WMO ' platchar])
    xlabel('Salinity');
    ylabel('Potential Temperature');
%AW    psname = ['ts_' platchar '.ps'];
%AW    eval(['print -dpsc ' psname])
  %}  

    
