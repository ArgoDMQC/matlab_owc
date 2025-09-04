
function [ pa_grid_sal, pa_grid_ptmp, pa_grid_pres, pa_grid_lat, pa_grid_long, pa_grid_dates , loaded_bbox ] = ...
    retr_region_ow( pa_wmo_numbers, pa_float_name, po_config_data, index, P, map_p_delta , loaded_bbox )

%
% function [ pa_grid_sal, pa_grid_ptmp, pa_grid_pres, pa_grid_lat, pa_grid_long, pa_grid_dates ] = retr_region_ow( pa_wmo_numbers, pa_float_name, po_config_data, index, P, map_p_delta ) ;
%
% This function retrieves historical data from the WMO boxes corresponding 
% to best hist observations, merges the CTD, BOT and Argo files, makes sure 
% longitude is continuous around the 0-360 degree mark, and converts the dates of the
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
%
% % Delphine Dobler (DD), August 2024: 
%     - 3.1: correct for iteration when there is only one box to load 
%           (consecutive to update in the indexation for performance 
%           improvement); adapt the creation of output vector to the
%           new index input.
%     - 3.2: save loaded data, and at each following profile, 
%            reload only new boxes and delete unused.



zz = (isnan(P)==0);
max_p = max(P(zz))+map_p_delta; % max depth to retrieve historical data is deepest Argo point plus map_p_delta

pa_index_0 = 0;

pa_grid_sal   = [ ] ;
pa_grid_ptmp  = [ ] ;
pa_grid_pres  = [ ] ;
pa_grid_lat   = [ ] ;
pa_grid_long  = [ ] ;
pa_grid_dates = [ ] ;

%[ max_depth, how_many_cols ] = size(pa_grid_pres);

% DD (2024/08-3.2): 
% Clean up loaded_bbox for unused boxes for the treated cycle (memory limits)
list_CTD_loaded_boxes=fields(loaded_bbox.CTD);
list_ARGO_loaded_boxes=fields(loaded_bbox.ARGO);
list_BOTTLE_loaded_boxes=fields(loaded_bbox.BOTTLE);

if size(list_CTD_loaded_boxes,1) > 1
    
    list_CTD_loaded_boxes=strrep(list_CTD_loaded_boxes,'box_','');
    list_CTD_loaded_boxes=str2double(list_CTD_loaded_boxes(2:end));
    list_CTD_loaded_boxes_to_delete=list_CTD_loaded_boxes(~ismember(list_CTD_loaded_boxes,pa_wmo_numbers(:,1)));
    
    if size(list_CTD_loaded_boxes_to_delete,1)> 0
        for i_box_to_delete=1:size(list_CTD_loaded_boxes_to_delete)
            loaded_bbox.CTD=rmfield(loaded_bbox.CTD,(['box_' num2str(list_CTD_loaded_boxes_to_delete(i_box_to_delete))]));
        end
    end
    
end

if size(list_ARGO_loaded_boxes,1) > 1
    
    list_ARGO_loaded_boxes=strrep(list_ARGO_loaded_boxes,'box_','');
    list_ARGO_loaded_boxes=str2double(list_ARGO_loaded_boxes(2:end));
    list_ARGO_loaded_boxes_to_delete=list_ARGO_loaded_boxes(~ismember(list_ARGO_loaded_boxes,pa_wmo_numbers(:,1)));
    
    if size(list_ARGO_loaded_boxes_to_delete)> 0
        for i_box_to_delete=1:size(list_ARGO_loaded_boxes_to_delete)
            loaded_bbox.ARGO=rmfield(loaded_bbox.ARGO,(['box_' num2str(list_ARGO_loaded_boxes_to_delete(i_box_to_delete))]));
        end
    end
    
end

if size(list_BOTTLE_loaded_boxes,1) > 1
    
    list_BOTTLE_loaded_boxes=strrep(list_BOTTLE_loaded_boxes,'box_','');
    list_BOTTLE_loaded_boxes=str2double(list_BOTTLE_loaded_boxes(2:end));
    list_BOTTLE_loaded_boxes_to_delete=list_BOTTLE_loaded_boxes(~ismember(list_BOTTLE_loaded_boxes,pa_wmo_numbers(:,1)));
    
    if size(list_BOTTLE_loaded_boxes_to_delete)> 0
        for i_box_to_delete=1:size(list_BOTTLE_loaded_boxes_to_delete)
            loaded_bbox.BOTTLE=rmfield(loaded_bbox.BOTTLE,(['box_' num2str(list_BOTTLE_loaded_boxes_to_delete(i_box_to_delete))]));
        end
    end
    
end


% % DD (2024/08-3.1): correct for loop iteration when there is only one box to load.
if size(pa_wmo_numbers,1) == 4
    n_boxes=size(pa_wmo_numbers,2);
else
    n_boxes=size(pa_wmo_numbers,1);
end
%for ln_index = 1:length(pa_wmo_numbers)
for ln_index = 1:n_boxes
    
    for ntyp = 2:4
        
        if( ~isnan(pa_wmo_numbers(ln_index,1)) & pa_wmo_numbers(ln_index,ntyp) ) % check to see if we are supposed to load this data type
            
            if ntyp == 2 % the 2nd column denotes CTD data
                
                % Test if it is already loaded % change DD (2024/08-3.2)
                A=isfield(loaded_bbox.CTD,['box_' num2str(pa_wmo_numbers(ln_index,1))]);
                if A == false
                    lo_box_data = load( strcat( po_config_data.HISTORICAL_DIRECTORY, po_config_data.HISTORICAL_CTD_PREFIX, sprintf( '%4d', pa_wmo_numbers(ln_index,1))));
                    loaded_bbox.CTD.(['box_' num2str(pa_wmo_numbers(ln_index,1))])=lo_box_data;
                else
                    lo_box_data=loaded_bbox.CTD.(['box_' num2str(pa_wmo_numbers(ln_index,1))]);
                end

            elseif ntyp == 3 % the 3rd column denotes historical data
                
                % Test if it is already loaded % change DD (2024/08-3.2)
                A=isfield(loaded_bbox.BOTTLE,['box_' num2str(pa_wmo_numbers(ln_index,1))]);
                if A == false
                    lo_box_data = load( strcat( po_config_data.HISTORICAL_DIRECTORY, po_config_data.HISTORICAL_BOTTLE_PREFIX, sprintf( '%4d', pa_wmo_numbers(ln_index,1))));
                     % In some cases the dimensions of the data may not match, causing issues with indexing below 
                    if numel(lo_box_data.lat)==numel(lo_box_data.pres) && size(lo_box_data.lat,2)~=size(lo_box_data.pres,2)
                        lo_box_data.pres = lo_box_data.pres';
                        lo_box_data.ptmp = lo_box_data.ptmp';
                        lo_box_data.sal = lo_box_data.sal';
                        lo_box_data.temp = lo_box_data.temp';
                    end
                    loaded_bbox.BOTTLE.(['box_' num2str(pa_wmo_numbers(ln_index,1))])=lo_box_data;
                else
                    lo_box_data=loaded_bbox.BOTTLE.(['box_' num2str(pa_wmo_numbers(ln_index,1))]);
                    
                end
            elseif ntyp == 4 % the 4th column denotes Argo data
                
                % Test if it is already loaded % change DD (2024/08-3.2)
                A=isfield(loaded_bbox.ARGO,['box_' num2str(pa_wmo_numbers(ln_index,1))]);
                if A == false
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
                
                    loaded_bbox.ARGO.(['box_' num2str(pa_wmo_numbers(ln_index,1))])=lo_box_data;
                else
                    lo_box_data=loaded_bbox.ARGO.(['box_' num2str(pa_wmo_numbers(ln_index,1))]);
                end
                
                %-----------------------------------------------------------------------

            end


            % DD (2024/08-3.1)
            % index is a concatenation of all stations inside the ellipse
            % whatever the initial box.
            % pa_index is a reconstruction of begin and end index of the
            % current box.
            i1 = length(lo_box_data.lat);
            i2 = 1:i1;
            pa_index = pa_index_0 + i2; % indices for this set of data
            pa_index_0 = pa_index_0 + i1; % increment beginning index
            res=ismember(pa_index,index);
            i_ok=find(res==1);
            n_clim_profiles_to_load=length(i_ok);
            %disp(['box_id=',num2str(pa_wmo_numbers(ln_index,1)),...
            %    ',ntyp=',num2str(ntyp),',n_clim_profiles_to_load=' num2str(n_clim_profiles_to_load)])
            
            if n_clim_profiles_to_load > 0 % not really necessary anymore, only for safe-keeping
                
                % no need to re-scan for each
                % station as the index is now saved. The extraction is also
                % treated with array operations, instead of looping.
                tmp_pa_grid_pres  = lo_box_data.pres(:,i_ok) ;
                tmp_pa_grid_sal   = lo_box_data.sal(:,i_ok)  ; 
                tmp_pa_grid_ptmp  = lo_box_data.ptmp(:,i_ok) ;
                
                too_deep=tmp_pa_grid_pres>max_p;
                tmp_pa_grid_pres(too_deep) = NaN;
                tmp_pa_grid_sal(too_deep)  = NaN;
                tmp_pa_grid_ptmp(too_deep) = NaN;
                
                i_entire_NaN_rows=all(isnan(tmp_pa_grid_pres),2);
                tmp_pa_grid_pres(i_entire_NaN_rows,:)=[];
                tmp_pa_grid_sal(i_entire_NaN_rows,:)=[];
                tmp_pa_grid_ptmp(i_entire_NaN_rows,:)=[];
                
                [n_pres_level_1,n_obs_1]=size(pa_grid_pres);
                [n_pres_level_2,n_obs_2]=size(tmp_pa_grid_pres);
                d_pres_level=n_pres_level_2-n_pres_level_1;
                if n_pres_level_1 > 0
                    if d_pres_level > 0
                        nan_array=NaN*ones(d_pres_level,n_obs_1);
                        pa_grid_pres=cat(2,cat(1,pa_grid_pres,nan_array), tmp_pa_grid_pres);
                        pa_grid_sal =cat(2,cat(1,pa_grid_sal ,nan_array), tmp_pa_grid_sal );
                        pa_grid_ptmp=cat(2,cat(1,pa_grid_ptmp,nan_array), tmp_pa_grid_ptmp);

                    else
                        if d_pres_level < 0
                            nan_array=NaN*ones(-d_pres_level,n_obs_2);
                            pa_grid_pres=cat(2,pa_grid_pres,cat(1,tmp_pa_grid_pres,nan_array));
                            pa_grid_sal =cat(2,pa_grid_sal ,cat(1,tmp_pa_grid_sal ,nan_array));
                            pa_grid_ptmp=cat(2,pa_grid_ptmp,cat(1,tmp_pa_grid_ptmp,nan_array));
                        else
                            pa_grid_pres=[pa_grid_pres, tmp_pa_grid_pres];
                            pa_grid_sal =[pa_grid_sal , tmp_pa_grid_sal ];
                            pa_grid_ptmp=[pa_grid_ptmp, tmp_pa_grid_ptmp];
                        end

                    end
                else
                    pa_grid_pres=[pa_grid_pres, tmp_pa_grid_pres];
                    pa_grid_sal =[pa_grid_sal , tmp_pa_grid_sal ];
                    pa_grid_ptmp=[pa_grid_ptmp, tmp_pa_grid_ptmp];
                end

                pa_grid_lat   = [pa_grid_lat   , lo_box_data.lat(i_ok)   ];
                pa_grid_long  = [pa_grid_long  , lo_box_data.long(i_ok)  ];
                pa_grid_dates = [pa_grid_dates , lo_box_data.dates(i_ok) ];    
                
            end %  if n_clim_profiles_to_load > 0
            % end of DD (2024/08-3.1)
                
%             %for n=1:i1  % load each station
%             
%                 %ok = find(index == pa_index(n));
%                 %if ~isempty(ok) % this entry is in the index list
% %aw                if n > length(lo_box_data.lat); keyboard; end
% 
%                   jj = find( isnan(lo_box_data.pres(:,n))==0 ); % only retrieve non-NaN values
%                   pres = lo_box_data.pres(jj,n);
%                   sal = lo_box_data.sal(jj,n);
%                   ptmp = lo_box_data.ptmp(jj,n);
% 
%                   too_deep = find(pres>max_p);
%                   pres(too_deep)=[];
%                   sal(too_deep)=[];
%                   ptmp(too_deep)=[];
% 
%                   new_depth = length(pres);
%                   if( new_depth>max_depth & max_depth~=0 ) % patch up number of rows
%                     pa_grid_pres = [ pa_grid_pres; ones( new_depth-max_depth, how_many_cols ) * NaN ];
%                     pa_grid_ptmp = [ pa_grid_ptmp; ones( new_depth-max_depth, how_many_cols ) * NaN ];
%                     pa_grid_sal = [ pa_grid_sal; ones( new_depth-max_depth, how_many_cols ) * NaN ];
%                   elseif( new_depth < max_depth )
%                     pres = [ pres; ones( max_depth-new_depth, 1 ) * NaN ];
%                     sal  = [ sal ; ones( max_depth-new_depth, 1 ) * NaN ];
%                     ptmp = [ ptmp; ones( max_depth-new_depth, 1 ) * NaN ];
%                   end
% 
%                   pa_grid_sal   = [ pa_grid_sal, sal ] ;
%                   pa_grid_ptmp  = [ pa_grid_ptmp, ptmp ] ;
%                   pa_grid_pres  = [ pa_grid_pres, pres ] ;
% 
%                   pa_grid_lat   = [ pa_grid_lat, lo_box_data.lat(n) ] ;
%                   pa_grid_long  = [ pa_grid_long, lo_box_data.long(n) ] ;
%                   pa_grid_dates = [ pa_grid_dates, lo_box_data.dates(n) ] ;
% 
%                   [ max_depth, how_many_cols ] = size(pa_grid_pres);

%                 end
%             end % for n=1:i1

        end % if( ~isnan(pa_wmo_numbers(ln_index,1)) & pa_wmo_numbers(ln_index,ntyp) )
    end %ntyp
end % for ln_index = 1:size(pa_wmo_numbers,1)


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

