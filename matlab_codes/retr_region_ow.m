
function [ pa_grid_sal, pa_grid_ptmp, pa_grid_pres, pa_grid_lat, pa_grid_long, pa_grid_dates ] = retr_region_ow( pa_wmo_numbers, pa_float_name, po_config_data, index, P, map_p_delta ) ;

%
% function [ pa_grid_sal, pa_grid_ptmp, pa_grid_pres, pa_grid_lat, pa_grid_long, pa_grid_dates ] = retr_region_ow( pa_wmo_numbers, pa_float_name, po_config_data, index, P, map_p_delta ) ;
%
% This function retrieves historical data from the 25 selected WMO boxes,
% merges the CTD, BOT and Argo files, makes sure longitude is continuous
% around the 0-360 degree mark, and converts the dates of the
% historical data from time format2 to year.mo+ (changedates.m).
%
% pa_wmo_numbers can be NaN: when float profiles are out of range (65N, 65S),
% or when there's no .mat file in that box, denoted by 0 (e.g. box on land).
%
% Historical data have lat,long and dates organised in single rows.
% The output from this function gives lat,long and dates in columns,
% just for ease of checking .... really doesn't matter.
%
% Annie Wong, December 2007.
%
% modified from get_region.m, Dec 2006, Breck Owens.
%
% Isabelle Gaboury, 26 Sep. 2017: Added check on the dimensions of bottle data.



zz = find(isnan(P)==0);
max_p = max(P(zz))+map_p_delta; % max depth to retrieve historical data is deepest Argo point plus map_p_delta

pa_index_0 = 0;

pa_grid_sal   = [ ] ;
pa_grid_ptmp  = [ ] ;
pa_grid_pres  = [ ] ;
pa_grid_lat   = [ ] ;
pa_grid_long  = [ ] ;
pa_grid_dates = [ ] ;

[ max_depth, how_many_cols ] = size(pa_grid_pres);


for ln_index = 1:length(pa_wmo_numbers)
    for ntyp = 2:4
        if( ~isnan(pa_wmo_numbers(ln_index,1)) & pa_wmo_numbers(ln_index,ntyp) ) % check to see if we are supposed to load this data type
            if ntyp == 2 % the 2nd column denotes CTD data
                lo_box_data = load( strcat( po_config_data.HISTORICAL_DIRECTORY, po_config_data.HISTORICAL_CTD_PREFIX, sprintf( '%4d', pa_wmo_numbers(ln_index,1))));
            elseif ntyp == 3 % the 3rd column denotes historical data
                lo_box_data = load( strcat( po_config_data.HISTORICAL_DIRECTORY, po_config_data.HISTORICAL_BOTTLE_PREFIX, sprintf( '%4d', pa_wmo_numbers(ln_index,1))));
                 % In some cases the dimensions of the data may not match, causing issues with indexing below 
                if numel(lo_box_data.lat)==numel(lo_box_data.pres) && size(lo_box_data.lat,2)~=size(lo_box_data.pres,2)
                    lo_box_data.pres = lo_box_data.pres';
                    lo_box_data.ptmp = lo_box_data.ptmp';
                    lo_box_data.sal = lo_box_data.sal';
                    lo_box_data.temp = lo_box_data.temp';
                end
            elseif ntyp == 4 % the 4th column denotes Argo data
                lo_box_data = load( strcat( po_config_data.HISTORICAL_DIRECTORY, po_config_data.HISTORICAL_ARGO_PREFIX, sprintf( '%4d', pa_wmo_numbers(ln_index,1))));
                % exclude Argo float being analysed from the Argo reference data selection,
                % must do this step before concatenating the vectors, because "index" comes
                % from "get_region_ow.m", which includes this step ---------------------
       	        not_use=[];
                for i=1:length(lo_box_data.lat)
                  profile=lo_box_data.source{i};
                  jj=findstr(profile,'_');
                  ref_float=profile(1:jj-1);
                  kk=findstr(pa_float_name, ref_float);
                  if(isempty(kk)==0)
                    not_use=[not_use,i];
                  end
                end
                lo_box_data.lat(not_use)=[];
                lo_box_data.long(not_use)=[];
                lo_box_data.dates(not_use)=[];
                lo_box_data.sal(:,not_use)=[];
                lo_box_data.ptmp(:,not_use)=[];
                lo_box_data.pres(:,not_use)=[];
                %-----------------------------------------------------------------------

            end

            % check index to see if this is a station that should be loaded
            i1 = length(lo_box_data.lat);
            i2 = 1:i1;
            pa_index = pa_index_0 + i2; % indices for this set of data
            pa_index_0 = pa_index_0 + i1; % increment beginning index

            for n=1:i1  % load each station
                ok = find(index == pa_index(n));
                if ~isempty(ok) % this entry is in the index list
%aw                if n > length(lo_box_data.lat); keyboard; end

 		  jj = find( isnan(lo_box_data.pres(:,n))==0 ); % only retrieve non-NaN values
                  pres = lo_box_data.pres(jj,n);
                  sal = lo_box_data.sal(jj,n);
                  ptmp = lo_box_data.ptmp(jj,n);

                  too_deep = find(pres>max_p);
                  pres(too_deep)=[];
                  sal(too_deep)=[];
                  ptmp(too_deep)=[];

                  new_depth = length(pres);
                  if( new_depth>max_depth & max_depth~=0 ) % patch up number of rows
                    pa_grid_pres = [ pa_grid_pres; ones( new_depth-max_depth, how_many_cols ) * NaN ];
                    pa_grid_ptmp = [ pa_grid_ptmp; ones( new_depth-max_depth, how_many_cols ) * NaN ];
                    pa_grid_sal = [ pa_grid_sal; ones( new_depth-max_depth, how_many_cols ) * NaN ];
                  elseif( new_depth < max_depth )
                    pres = [ pres; ones( max_depth-new_depth, 1 ) * NaN ];
                    sal  = [ sal ; ones( max_depth-new_depth, 1 ) * NaN ];
                    ptmp = [ ptmp; ones( max_depth-new_depth, 1 ) * NaN ];
                  end

                  pa_grid_sal   = [ pa_grid_sal, sal ] ;
                  pa_grid_ptmp  = [ pa_grid_ptmp, ptmp ] ;
                  pa_grid_pres  = [ pa_grid_pres, pres ] ;

                  pa_grid_lat   = [ pa_grid_lat, lo_box_data.lat(n) ] ;
                  pa_grid_long  = [ pa_grid_long, lo_box_data.long(n) ] ;
                  pa_grid_dates = [ pa_grid_dates, lo_box_data.dates(n) ] ;

                  [ max_depth, how_many_cols ] = size(pa_grid_pres);

                end
	    end %n=1:i1
        end
    end %ntyp
end %ln_index


% longitude goes from 0 to 360 degrees

ln_jj = find( pa_grid_long < 0 ) ;
pa_grid_long( ln_jj ) = 360 + pa_grid_long( ln_jj ) ;


% make sure longitude is continuous around the 0-360 degree mark

ln_kk = find( pa_grid_long>=320 & pa_grid_long<=360 ) ;
if( isempty( ln_kk ) == 0 )
   ln_ll = find( pa_grid_long>=0 & pa_grid_long<=40 ) ;
   pa_grid_long( ln_ll ) = 360 + pa_grid_long( ln_ll ) ;
end


% make pa_grid_sal, pa_grid_ptmp, pa_grid_pres have the same NaNs

ln_ii = find( isnan( pa_grid_sal ) == 1 ) ;
pa_grid_pres( ln_ii ) = NaN.*ones( 1, length( ln_ii ) ) ;
pa_grid_ptmp( ln_ii ) = NaN.*ones( 1, length( ln_ii ) ) ;
ln_ii = find( isnan( pa_grid_pres ) == 1 ) ;
pa_grid_sal( ln_ii ) = NaN.*ones( 1, length( ln_ii ) ) ;
pa_grid_ptmp( ln_ii ) = NaN.*ones( 1, length( ln_ii ) ) ;
ln_ii = find( isnan( pa_grid_ptmp ) == 1 ) ;
pa_grid_sal( ln_ii ) = NaN.*ones( 1, length( ln_ii ) ) ;
pa_grid_pres( ln_ii ) = NaN.*ones( 1, length( ln_ii ) ) ;

pa_grid_dates = changedates( pa_grid_dates ) ;


% turns rows into columns

pa_grid_lat = pa_grid_lat' ;
pa_grid_long = pa_grid_long' ;

