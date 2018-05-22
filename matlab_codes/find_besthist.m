
function [ index ] = find_besthist( grid_lat, grid_long, grid_dates, grid_Z, LAT, LONG, DATES, Z, latitude_large, latitude_small, longitude_large, longitude_small, phi_large, phi_small, age, age_large, map_pv_use, max_casts )

%
% Find ln_max_casts unique historical points that are most strongly correlated with the float profile.
%
% Rewritten Dec 2006 to only use lat, lon, date and water depth as input
% and to return index into station list.
%
% Annie Wong, June 2010.
% Breck Owens, December 2006.
%
% Cecile Cabanes, June 2013
% track "change config 129" : add age_large when computing correlation_large



% avoid having PV=0, because if PV=0, ellipse>1, and no historical points will be selected ---

PV_float = (2*7.292*10^-5.*sin(LAT.*pi/180))./Z;
PV_hist = (2*7.292*10^-5.*sin(grid_lat.*pi/180))./grid_Z;

if(PV_float==0)PV_float=1*10^-5;end
jj=find(PV_hist==0);
PV_hist(jj)=1*10^-5;


% pick out points within the ellipse ---

if( map_pv_use==1 ) % if PV is wanted
  ellipse = sqrt( (grid_long-LONG).^2./(longitude_large*3).^2 + (grid_lat-LAT).^2./(latitude_large*3).^2 +...
     ((PV_float-PV_hist)./sqrt( PV_float.^2+PV_hist.^2 )./phi_large).^2 ) ;
else % if PV is unwanted ---
  ellipse = sqrt((grid_long-LONG).^2./(longitude_large*3).^2 + (grid_lat-LAT).^2./(latitude_large*3).^2) ;
end

index_ellipse = find( ellipse<1 ) ;
hist_long  = grid_long (index_ellipse) ;
hist_lat   = grid_lat  (index_ellipse) ;
hist_dates = grid_dates(index_ellipse) ;
hist_Z = grid_Z(index_ellipse) ;

index=index_ellipse; %if <max_casts


if(length(index_ellipse)>max_casts)

%% pick max_casts/3 random points ---
    rng(length(index_ellipse)); % ensures that the result will be the same each time the mapping is done with the same reference dataset

    index_random = round(rand(1,ceil(max_casts/3)).*length(hist_long));
    kk=find(index_random==0);
    index_random(kk)=ones(1,length(kk));
    index_random=unique(index_random);

    index_remain = [1:length(hist_long)];
    ii=[];
    for h=1:length(hist_long)
        a=find(index_random==index_remain(h));
        if(isempty(a)==0)
            ii=[ii,h];
        end
    end
    index_remain(ii)=[];


%% sort remaining points by large spatial correlations ---

    remain_hist_lat = hist_lat(index_remain) ;
    remain_hist_long = hist_long(index_remain) ;
    remain_hist_dates = hist_dates(index_remain) ;
    remain_hist_Z = hist_Z(index_remain) ;

    PV_hist = (2*7.292*10^-5.*sin(remain_hist_lat.*pi/180))./remain_hist_Z;

    if( map_pv_use==1 ) % if PV is wanted
        %correlation_large = (remain_hist_long-LONG).^2./longitude_large.^2 + (remain_hist_lat-LAT).^2./latitude_large.^2 +...
        %    ( (PV_float-PV_hist)./sqrt( PV_float.^2+PV_hist.^2 )./phi_large ).^2 ;
        correlation_large = (remain_hist_long-LONG).^2./longitude_large.^2 + (remain_hist_lat-LAT).^2./latitude_large.^2 +...
            (remain_hist_dates-DATES).^2./age_large.^2 + ( (PV_float-PV_hist)./sqrt( PV_float.^2+PV_hist.^2 )./phi_large ).^2 ;
        
    else % if PV is unwanted ---
       % correlation_large = (remain_hist_long-LONG).^2./longitude_large.^2 + (remain_hist_lat-LAT).^2./latitude_large.^2 ;
        correlation_large = (remain_hist_long-LONG).^2./longitude_large.^2 + (remain_hist_lat-LAT).^2./latitude_large.^2 +...
            (remain_hist_dates-DATES).^2./age_large.^2 ;
    end

    [ sorted_correlation_large, index_large ] = sort( correlation_large ) ;

    remain_hist_lat = remain_hist_lat( index_large ) ;
    remain_hist_long = remain_hist_long( index_large ) ;
    remain_hist_dates = remain_hist_dates( index_large ) ;
    remain_hist_Z = remain_hist_Z( index_large ) ;


%% sort remaining points by short spatial and temporal correlations ---

    lsegment2 = ceil(max_casts/3) + ceil(max_casts/3) - length(index_random); %length of segment2 is 1/3 of max_casts plus whatever is deficient from unique(random)

    index_z = lsegment2+1 : length(remain_hist_long);
    zlong  = remain_hist_long ( index_z ) ;
    zlat   = remain_hist_lat  ( index_z ) ;
    zdates = remain_hist_dates( index_z ) ;
    zZ = remain_hist_Z( index_z ) ;

    PV_hist = (2*7.292*10^-5.*sin(zlat.*pi/180))./zZ;

    if( map_pv_use==1 ) % if PV is wanted
        correlation_small = (zlong-LONG).^2./longitude_small.^2 + (zlat-LAT).^2./latitude_small.^2 +...
            (zdates-DATES).^2./age.^2 + ( (PV_float-PV_hist)./sqrt( PV_float.^2+PV_hist.^2 )./phi_small ).^2 ;
    else % if PV is unwanted
        correlation_small = (zlong-LONG).^2./longitude_small.^2 + (zlat-LAT).^2./latitude_small.^2 +...
            (zdates-DATES).^2./age.^2 ;
    end

    [ sorted_correlation_small, index_small ] = sort( correlation_small ) ;


%% piece the 3 steps together ---

    leftover = max_casts - length(index_random) - lsegment2; %same as leftover=max_casts-2*ceil(max_casts/3);

    % get index into original list of stations
    i1 = index_ellipse( index_random );
    % index for hist_lat(index_random)
    i2 = index_ellipse( index_remain( index_large( 1:lsegment2 ) ) ); % index for remain_hist_lat(1:lsegment2)
    i3 = index_ellipse( index_remain( index_large( index_z( index_small( 1:leftover ) ) ) ) ); % index for zlat
    index = [i1; i2; i3];
    index = sort(index);

end %if(length(index_ellipse)>max_casts)


