
function update_salinity_mapping( pn_float_dir, pn_float_name, po_system_configuration )

% function update_salinity_mapping( pn_float_dir, pn_float_name, po_system_configuration )
%
% Annie Wong, June 2010.
% Breck Owens, November 2007.

% pn_float_dir='testfloats/';
% pn_float_name='R6900395';
% po_system_configuration = load_configuration( 'ow_config.txt' );


% Cecile Cabanes, June 2013 : use of "map_age_large" the time scale used to mapped the large scale field
% to track modifications: see "changes config 129"  in this program
%
% Giulio Notarstefano, Dec 2018: modify the threshold values for the
% quality control of mapped salinities in the Black Sea
%
% Annie Wong, June 2020: rename "scale_age" to "scale_age_small", and create "scale_age_large" with NaN for old OW mapped profiles.
% This is to make the temporal scale variable names consistent with the lat/long scale variable names in OWC.
%
% Delphine Dobler (DD), August 2024: 
%            1 - Auto-create la_wmo_boxes from climatological repositores when la_wmo_boxes does not exist
%

% load float source data ----------------------------------------

filename = fullfile( po_system_configuration.FLOAT_SOURCE_DIRECTORY, pn_float_dir, strcat( pn_float_name, po_system_configuration.FLOAT_SOURCE_POSTFIX ) );
lo_float_source_data = load( filename );

PROFILE_NO = lo_float_source_data.PROFILE_NO ;
[ ln_float_level_count, ln_float_profile_count ] = size( lo_float_source_data.SAL ) ;


% load mapping parameters -------
% DD : 2024/08-1
la_wmo_boxes_file=strcat( po_system_configuration.CONFIG_DIRECTORY, po_system_configuration.CONFIG_WMO_BOXES);
if ~exist(la_wmo_boxes_file)
    create_la_wmo_boxes_file(po_system_configuration)
end
load( la_wmo_boxes_file, 'la_wmo_boxes' ) ;

ln_max_casts = str2num( po_system_configuration.CONFIG_MAX_CASTS ) ;
map_use_pv = str2num( po_system_configuration.MAP_USE_PV );
map_use_saf = str2num( po_system_configuration.MAP_USE_SAF );
longitude_large = str2double( po_system_configuration.MAPSCALE_LONGITUDE_LARGE);
longitude_small = str2double( po_system_configuration.MAPSCALE_LONGITUDE_SMALL);
latitude_large  = str2double( po_system_configuration.MAPSCALE_LATITUDE_LARGE );
latitude_small  = str2double( po_system_configuration.MAPSCALE_LATITUDE_SMALL );
phi_large = str2double( po_system_configuration.MAPSCALE_PHI_LARGE );
phi_small = str2double( po_system_configuration.MAPSCALE_PHI_SMALL );
map_age_small = str2double( po_system_configuration.MAPSCALE_AGE_SMALL );  %AW June 2020
%% added CC
map_age_large = str2double( po_system_configuration.MAPSCALE_AGE_LARGE );
%%
map_p_delta = str2double( po_system_configuration.MAP_P_DELTA );
map_p_exclude = str2double( po_system_configuration.MAP_P_EXCLUDE );


% load precalculated mapped data --------

ls_float_mapped_filename = fullfile( po_system_configuration.FLOAT_MAPPED_DIRECTORY, pn_float_dir, strcat( po_system_configuration.FLOAT_MAPPED_PREFIX, pn_float_name, po_system_configuration.FLOAT_MAPPED_POSTFIX ) );

if exist( ls_float_mapped_filename, 'file' )
   load( ls_float_mapped_filename );
   ln_profile_index = size( la_mapped_sal, 2 );
   
   % check to see if this an older run without saf
   if ~exist('use_saf', 'var')
      use_saf = zeros(size(use_pv));
   end

   % check to see if this is an older run with scale_age and no scale_age_large --- AW June 2020
   if exist('scale_age', 'var')
      scale_age_small = scale_age;
      clear scale_age
   end

   if ~exist('scale_age_large', 'var')
      scale_age_large = NaN.*ones(1,length(scale_age_small));
   end
   % check to see if this is an older run with scale_age and no scale_age_large --- AW June 2020

   [ max_depth, how_many_cols ] = size(la_mapped_sal);
   new_depth = ln_float_level_count;
   if( new_depth>max_depth & max_depth~=0 ) % patch up number of rows
     la_mapped_sal = [ la_mapped_sal; ones( new_depth-max_depth, how_many_cols ) * NaN ];
     la_mapsalerrors = [ la_mapsalerrors; ones( new_depth-max_depth, how_many_cols ) * NaN ];
     la_noise_sal = [ la_noise_sal; ones( new_depth-max_depth, how_many_cols ) * NaN ];
     la_signal_sal = [ la_signal_sal; ones( new_depth-max_depth, how_many_cols ) * NaN ];
     la_ptmp = [ la_ptmp; ones( new_depth-max_depth, how_many_cols ) * NaN ];
   end
else
   ln_profile_index = 0 ;
   la_profile_no = NaN ;
   selected_hist = [] ;
end


% compare profile numbers in float source matrix and mapped data matrix ------

missing_profile_index = [] ;

for i = 1 : length( PROFILE_NO )
  a = find( la_profile_no==PROFILE_NO(i) ) ;
  if( isempty(a)==1 )
     missing_profile_index = [ missing_profile_index, i ] ;
  end
end
clear a i


% update mapped data matrix by missing_profile_index --------------------

for i = 1 : length( missing_profile_index )
  tic
  disp(['UPDATE_SALINITY_MAPPING: Working on profile ' num2str(i)])

  j = missing_profile_index( i ) ;
  ln_profile_index = ln_profile_index + 1 ;

  % append profile numbers

  la_profile_no( ln_profile_index ) = PROFILE_NO( j ) ;

  % initialize output variables

  la_ptmp ( :, ln_profile_index ) = NaN.*ones( ln_float_level_count, 1 ) ;
  la_mapped_sal  ( :, ln_profile_index ) = NaN.*ones( ln_float_level_count, 1 ) ;
  la_mapsalerrors( :, ln_profile_index ) = NaN.*ones( ln_float_level_count, 1 ) ;
  la_noise_sal   ( :, ln_profile_index ) = NaN.*ones( ln_float_level_count, 1 ) ;
  la_signal_sal  ( :, ln_profile_index ) = NaN.*ones( ln_float_level_count, 1 ) ;
  scale_long_large( ln_profile_index ) = NaN;
  scale_lat_large ( ln_profile_index ) = NaN;
  scale_long_small( ln_profile_index ) = NaN;
  scale_lat_small ( ln_profile_index ) = NaN;
  scale_phi_large ( ln_profile_index ) = NaN;
  scale_phi_small ( ln_profile_index ) = NaN;
  scale_age_small ( ln_profile_index ) = NaN; %AW June 2020
  scale_age_large ( ln_profile_index ) = NaN;
  use_pv ( ln_profile_index ) = NaN;
  use_saf ( ln_profile_index ) = NaN;
  p_delta ( ln_profile_index ) = NaN;
  p_exclude ( ln_profile_index ) = NaN;

  % get float data

  LAT = lo_float_source_data.LAT( j ) ;
  LONG = lo_float_source_data.LONG( j ) ;
  DATES = lo_float_source_data.DATES(j ) ;
  SAL = lo_float_source_data.SAL(:,j) ;
  PTMP = lo_float_source_data.PTMP(:,j) ;
  PRES = lo_float_source_data.PRES(:,j) ;
  
    if( isnan(LONG)==0 & isnan(LAT)==0 & isempty(find(isnan(PRES)==0))~=1 )

    if(LONG>180) % m_tbase inputs longitudes from 0 to +/- 180
       LONG1=LONG-360;
    else
       LONG1=LONG;
    end
    m_proj('mercator','long', [min(LONG1)-1, max(LONG1)+1], 'lat', [min(LAT)-1, max(LAT)+1] );
    [elev,x,y] = m_tbase( [min(LONG1)-1, max(LONG1)+1, min(LAT)-1, max(LAT)+1] );
    Z = -interp2( x,y,elev, LONG1, LAT, 'linear'); % -ve bathy values

    % get historical data

    [ la_wmo_numbers ] = find_25boxes( LONG, LAT, la_wmo_boxes ) ;
    % load in the positions and dates for the surrounding wmo boxes
    [ la_grid_lat, la_grid_long, la_grid_dates ] = get_region_ow( la_wmo_numbers, pn_float_name, po_system_configuration ) ;

    if( la_grid_lat~=999 ) % if no historical data is assigned to the float profile

        % calculate water depth for historical casts

        la_grid_long1=la_grid_long;
        gg=find(la_grid_long>180);
        la_grid_long1(gg)=la_grid_long(gg)-360; % m_tbase inputs longitudes from 0 to +/- 180

	m_proj('mercator','long', [min(la_grid_long1)-1, max(la_grid_long1)+1], 'lat', [min(la_grid_lat)-1, max(la_grid_lat)+1] );
        [elev,x,y] = m_tbase( [min(la_grid_long1)-1, max(la_grid_long1)+1, min(la_grid_lat)-1, max(la_grid_lat)+1] );
        la_grid_Z = -interp2( x,y,elev, la_grid_long1, la_grid_lat, 'linear'); % -ve bathy values

        % make LONG compatiable with la_grid_long at the 0-360 mark

        LONG2=LONG;
	if(isempty(find(la_grid_long>360))==0)
           if(LONG>=0&LONG<=20)
             LONG2=LONG+360;
           end
        end

        % find ln_max_casts historical points that are most strongly correlated with the float profile
        
        %[ index ] = find_besthist( la_grid_lat, la_grid_long, la_grid_dates, la_grid_Z, LAT, LONG2, DATES, Z, latitude_large, latitude_small, longitude_large, longitude_small, phi_large, phi_small, map_age, map_use_pv, ln_max_casts );      % change config 129 
	[ index ] = find_besthist( la_grid_lat, la_grid_long, la_grid_dates, la_grid_Z, LAT, LONG2, DATES, Z, latitude_large, latitude_small, longitude_large, longitude_small, phi_large, phi_small, map_age_small, map_age_large, map_use_pv, ln_max_casts );  %AW June 2020
	
        clear la_grid_lat la_grid_long la_grid_dates
        
	[ la_bhist_sal, la_bhist_ptmp, la_bhist_pres, la_bhist_lat, la_bhist_long, la_bhist_dates ] = retr_region_ow( la_wmo_numbers, pn_float_name, po_system_configuration, index, PRES, map_p_delta ) ;
        la_bhist_Z = la_grid_Z(index);

        % include JB's SAF frontal separation criteria if map_use_saf==1

        if(map_use_saf==1)

          [ la_bhist_sal2, la_bhist_ptmp2, la_bhist_pres2, la_bhist_lat2, la_bhist_long2, la_bhist_dates2, la_bhist_Z2 ] = ...
           frontalConstraintSAF(la_bhist_sal, la_bhist_ptmp, la_bhist_pres, la_bhist_lat, la_bhist_long, la_bhist_dates, la_bhist_Z, LAT, LONG, PRES, PTMP, SAL, po_system_configuration);

          [a2 b2]= size(la_bhist_sal2);
          if (b2>5) % Use frontal separation only if there are at least 5 profiles
            la_bhist_sal=la_bhist_sal2;
            la_bhist_ptmp=la_bhist_ptmp2;
            la_bhist_pres=la_bhist_pres2;
            la_bhist_lat=la_bhist_lat2;
            la_bhist_long=la_bhist_long2;
            la_bhist_dates=la_bhist_dates2;
            la_bhist_Z=la_bhist_Z2;
          end

          clear la_bhist_sal2 la_bhist_ptmp2 la_bhist_pres2 la_bhist_lat2 la_bhist_long2 la_bhist_dates2 la_bhist_Z2

        end %if(map_use_saf==1)

        % make LONG compatiable with la_bhist_long at the 0-360 mark

        if(isempty(find(la_bhist_long>360))==0)
           if(LONG>=0&LONG<=20)
              LONG=LONG+360;
           end
        end

        % interpolate historical casts onto the float theta surfaces

        la_hist_interpsal=NaN.*ones( ln_float_level_count, length(la_bhist_lat) ); %length(la_bhist_lat) can be smaller than ln_max_casts
        la_hist_interppres=NaN.*ones( ln_float_level_count, length(la_bhist_lat) ); %length(la_bhist_lat) can be smaller than ln_max_casts

        [ la_hist_interpsal, la_hist_interppres ] = interp_climatology( la_bhist_sal, la_bhist_ptmp, la_bhist_pres, SAL, PTMP, PRES );
        clear la_bhist_sal la_bhist_ptmp la_bhist_pres la_wmo_numbers index

        % map one float theta level at a time below map_p_exclude

	for ln_level = 1:ln_float_level_count

           if( isnan(SAL(ln_level))==0 & PRES(ln_level)>=map_p_exclude )

	    % for each float theta level, only use the la_hist_interpsal that are not NaNs

            ln_max_hist_casts = find( isnan(la_hist_interpsal(ln_level,:))==0 ) ;
            la_hist_sal   = la_hist_interpsal  ( ln_level, ln_max_hist_casts ) ;
            la_hist_pres  = la_hist_interppres ( ln_level, ln_max_hist_casts ) ;
            la_hist_long  = la_bhist_long ( ln_max_hist_casts ) ;
            la_hist_lat   = la_bhist_lat  ( ln_max_hist_casts ) ;
            la_hist_dates = la_bhist_dates( ln_max_hist_casts ) ;
            la_hist_Z = la_bhist_Z( ln_max_hist_casts ) ;

            % pick out points within +/- map_p_delta dbar of float pressure at that ptlevel

            compare_pres = la_hist_pres - PRES(ln_level);
            ii=[];
            ii=find( abs(compare_pres) < map_p_delta );
            la_hist_sal   = la_hist_sal  (ii) ;
            la_hist_long  = la_hist_long (ii) ;
            la_hist_lat   = la_hist_lat  (ii) ;
            la_hist_dates = la_hist_dates(ii) ;
            la_hist_Z = la_hist_Z(ii) ;
            
            % map historical data to float profiles -------------------------
            % check for outliers
            
            mean_sal = mean(la_hist_sal);
            signal_sal = signal( la_hist_sal ) ;
            
            if signal_sal>0
                bad = find(abs(la_hist_sal-mean_sal)/sqrt(signal_sal) >3);
                
                if ~isempty(bad)
                    la_hist_sal(bad)   = [];
                    la_hist_long(bad)  = [];
                    la_hist_lat(bad)   = [];
                    la_hist_dates(bad) = [];
                    la_hist_Z(bad)     = [];
                end
            end
            
            
            % require more than one historical salinity profile with unique position for noise.m; and
            % require more than five historical salinity profiles to produce realistic mapping error estimates
            if ( (length(unique(la_hist_long)) > 1 | length(unique(la_hist_lat)) > 1) & length(la_hist_sal) > 5 )
                
                
                % use large length scales to map original data
                % pass NaN for map_age to exclude temporal covariance for large scale mapping
                
                noise_sal  = noise( la_hist_sal, la_hist_lat, la_hist_long ) ;
                signal_sal = signal( la_hist_sal ) ;
               
               % change config 129
               
%                  [a1,b1,c1,d1]...
%                  = map_data_grid( la_hist_sal, [ LAT, LONG, DATES, Z ], ...
%                           [ la_hist_lat, la_hist_long, la_hist_dates, la_hist_Z ], ...
%                           longitude_large, latitude_large, NaN, signal_sal, noise_sal, phi_large, map_use_pv ) ;


               [a1,b1,c1,d1]...
		= map_data_grid( la_hist_sal, [ LAT, LONG, DATES, Z ], ...
			 [ la_hist_lat, la_hist_long, la_hist_dates, la_hist_Z ], ...
			 longitude_large, latitude_large, map_age_large, signal_sal, noise_sal, phi_large, map_use_pv ) ;

               % use short length scales and temporal scales to map residuals

               la_residualsal1 = la_hist_sal - c1' ;
               la_signalresidualsal = signal( la_residualsal1 ) ;

               [a2,b2,c2,d2]...
                = map_data_grid( la_residualsal1, [ LAT, LONG, DATES, Z ], ...
                         [ la_hist_lat, la_hist_long, la_hist_dates, la_hist_Z ], ...
				 longitude_small, latitude_small, map_age_small, la_signalresidualsal, noise_sal, phi_small, map_use_pv ) ;  %AW June 2020

               la_ptmp( ln_level, ln_profile_index ) = PTMP(ln_level);
               la_mapped_sal(  ln_level, ln_profile_index ) = a1' + a2' ;
               %la_mapsalerrors(ln_level, ln_profile_index ) = b2' ; 
               la_mapsalerrors(ln_level, ln_profile_index ) = sqrt(b1'.*b1' + b2'*b2') ; %change config 129 ;
               la_noise_sal(   ln_level, ln_profile_index ) = noise_sal ;
               la_signal_sal(  ln_level, ln_profile_index ) = signal_sal ;
               scale_long_large( ln_profile_index ) = longitude_large ;
               scale_lat_large ( ln_profile_index ) = latitude_large ;
               scale_long_small( ln_profile_index ) = longitude_small ;
               scale_lat_small ( ln_profile_index ) = latitude_small ;
               scale_phi_large ( ln_profile_index ) = phi_large ;
               scale_phi_small ( ln_profile_index ) = phi_small ;
               scale_age_small ( ln_profile_index ) = map_age_small ;  %AW June 2020
               scale_age_large ( ln_profile_index ) = map_age_large ;
               use_pv ( ln_profile_index ) = map_use_pv ;
               use_saf ( ln_profile_index ) = map_use_saf ;
               p_delta ( ln_profile_index ) = map_p_delta ;
               p_exclude ( ln_profile_index ) = map_p_exclude ;

               % only save selected historical points to conserve computer space

               if(isempty(selected_hist)==1)
                  selected_hist=[la_hist_long(1),la_hist_lat(1),la_profile_no(ln_profile_index)];
               end
               for k = 1 : length(la_hist_long)
                  [m,n] = size(selected_hist);
                  b = [ la_hist_long(k), la_hist_lat(k) ];
                  c = selected_hist(:,1:2) - ones(m,1)*b;
                  d = find( abs(c(:,1))<1/60&abs(c(:,2))<1/60 ); %within 1 min, do not save
		  if( isempty(d)==1 )
                    selected_hist = [ selected_hist; [ la_hist_long(k), la_hist_lat(k), la_profile_no(ln_profile_index) ] ];
                  end
               end

            end %if( length(la_hist_sal)<=1 )
          end %if profile levels is within ptlevels, map
        end %for ln_level = 1:ln_number_used_levels
     end %if( la_grid_lat==999 )
  end %if(isnan(LONG)==0&isnan(LAT)==0)
  toc
  fclose('all');
end %for i = 1 : length( missing_profile_index )


% quality control - subst all mapped_sal < 30 and > 40 with NaNs ----------

ii = find( la_mapped_sal<30|la_mapped_sal>40 ) ;
if min(lo_float_source_data.LAT)>41 && max(lo_float_source_data.LAT)<48 && min(lo_float_source_data.LONG)>27 && max(lo_float_source_data.LONG)<42 
ii = find( la_mapped_sal<15|la_mapped_sal>25 ) ;  % GN -> Black Sea
end

la_mapped_sal(ii) = NaN.*ones( 1,length(ii) ) ;

% sort the mapped data matrix by profile numbers ------------

[y,ii] = sort( la_profile_no ) ;

la_ptmp = la_ptmp (:,ii) ;
la_mapped_sal   = la_mapped_sal  (:,ii) ;
la_mapsalerrors = la_mapsalerrors(:,ii) ;
la_noise_sal    = la_noise_sal   (:,ii) ;
la_signal_sal   = la_signal_sal  (:,ii) ;
scale_long_large= scale_long_large(ii) ;
scale_lat_large = scale_lat_large (ii) ;
scale_long_small= scale_long_small(ii) ;
scale_lat_small = scale_lat_small (ii) ;
scale_phi_large = scale_phi_large(ii) ;
scale_phi_small = scale_phi_small(ii) ;
scale_age_small = scale_age_small(ii) ;  %AW June 2020
scale_age_large = scale_age_large(ii) ;
use_pv = use_pv(ii) ;
use_saf = use_saf(ii) ;
p_delta = p_delta(ii) ;
p_exclude = p_exclude(ii) ;
la_profile_no = la_profile_no(ii) ;

if(isempty(selected_hist)==0)
  [y,ii] = sort( selected_hist(:,3) ) ;
  selected_hist = selected_hist(ii,:) ;
end

%keyboard
% save the relevant data ----------------

save( ls_float_mapped_filename, 'la_ptmp', 'la_mapped_sal', 'la_mapsalerrors', ...
      'scale_long_large', 'scale_lat_large', 'scale_long_small', 'scale_lat_small', ...
      'scale_phi_large', 'scale_phi_small', 'scale_age_small', 'scale_age_large', 'use_pv', 'use_saf', 'p_delta', 'p_exclude', ...
      'la_noise_sal', 'la_signal_sal', 'la_profile_no', 'selected_hist') ;  %AW June 2020


% clear memory ----

fclose('all');

