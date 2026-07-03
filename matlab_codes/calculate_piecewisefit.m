
function calculate_piecewisefit( pn_float_dir, pn_float_name, po_system_configuration )

%
% Annie Wong, October 2008
% Breck Owens, November 2007
% Cecile Cabanes, June 2013 : calculate off-diagonal terms for error estimate: add horizontal covariance   to track changes: see "change config 129"
% Delphine Dobler (DD), September 2024 : 
%            4 - save computed theta levels in a file 
% Delphine Dobler (DD), June 2026:
%     - correct bug in longitude domain change
%     - separate continuous piecewisefit segments and theta sets segments
%     - change breaks handling with piecewisefit segments: the code now
%     handles to which piecewisefit segment prescribed breaks belong to
%     - review check, info, warning and error messages for breaks and max_breaks
%     - remove unused parameters and systematicaaly allocate size for array
%     parameters

%pn_float_dir='testfloats/';
%pn_float_name='robbins4900178';
%po_system_configuration = load_configuration( 'ow_config.txt' );


% load data from /float_source and /float_mapped --------------

lo_float_source_data = load( strcat( po_system_configuration.FLOAT_SOURCE_DIRECTORY, pn_float_dir, pn_float_name, po_system_configuration.FLOAT_SOURCE_POSTFIX ) );

LAT = lo_float_source_data.LAT; % positions of the floats
LONG = lo_float_source_data.LONG;
DATES = lo_float_source_data.DATES;
SAL = lo_float_source_data.SAL; % salinity from the floats
PTMP = lo_float_source_data.PTMP; % potential temperature from the floats
PRES = lo_float_source_data.PRES; % pressure from the floats
PROFILE_NO = lo_float_source_data.PROFILE_NO; % profile number
x_in = repmat( PROFILE_NO, 10, 1);

lo_float_mapped_data = load( strcat( po_system_configuration.FLOAT_MAPPED_DIRECTORY, pn_float_dir, po_system_configuration.FLOAT_MAPPED_PREFIX, pn_float_name, po_system_configuration.FLOAT_MAPPED_POSTFIX ) );

mapped_sal = lo_float_mapped_data.la_mapped_sal; % salinity from climatology mapped to float locations and times
mapsalerrors = lo_float_mapped_data.la_mapsalerrors; % mapping errors for salinity
la_ptmp = lo_float_mapped_data.la_ptmp; % float potential temperature where mapping is done



%+++++ change config 129----FROM HERE

% retrieve coordinate (XYTZ) of the float position (coord_float) that is used in build_cov.m
%if ~isempty(lo_float_mapped_data.selected_hist)
% the following lines were added to make sure LONG, LAT, DATES are stored in
    % column vectors (from Dirk Slawinski-CSIRO) 
if size(LONG,1)>1 
	LONG=LONG';
end
if size(LAT,1)>1  
	LAT=LAT';
end
if size(DATES,1)>1
	DATES=DATES';
end

% Calculate elevation at the float position 
% DD 2026/06/30 - fixing longitude domain formula
LONG1=LONG;
LONG1(LONG1>180)=LONG1(LONG1>180)-360;
% end of fix

m_proj('mercator','long', [min(LONG1)-1, max(LONG1)+1], 'lat', [min(LAT)-1, max(LAT)+1] );
[elev,x,y] = m_tbase( [min(LONG1)-1, max(LONG1)+1, min(LAT)-1, max(LAT)+1] );
Z = -interp2( x,y,elev, LONG1, LAT, 'linear'); % -ve bathy values

coord_float = [LONG',LAT',DATES',Z'];

if isempty(lo_float_mapped_data.selected_hist)    
   warning('No reference data selected for any profiles')
end
%+++++ change config 129-----TO HERE


% load calibration settings -----------------

lo_float_calseries = load( strcat( po_system_configuration.FLOAT_CALIB_DIRECTORY, pn_float_dir, po_system_configuration.FLOAT_CALSERIES_PREFIX , pn_float_name, po_system_configuration.FLOAT_MAPPED_POSTFIX  ) );

calseries = lo_float_calseries.calseries; %for theta levels selection
calseries_for_cont_pwf = lo_float_calseries.calseries_for_cont_pwf; % for timeseries interpolation (piecewise linear fit)
max_breaks = lo_float_calseries.max_breaks;
breaks = lo_float_calseries.breaks;
use_theta_gt = lo_float_calseries.use_theta_gt;
use_theta_lt = lo_float_calseries.use_theta_lt;
use_pres_gt = lo_float_calseries.use_pres_gt;
use_pres_lt = lo_float_calseries.use_pres_lt;
use_percent_gt = lo_float_calseries.use_percent_gt;


%---- Breck's add-on -------

[m,n] = size(PRES);

cal_SAL = NaN*ones(m,n);
cal_SAL1 = NaN*ones(m,n);
cal_SAL_err = NaN*ones(m,n);
cal_COND = NaN*ones(m,n);
cal_COND_err = NaN*ones(m,n);
xfit = NaN*ones(1,n);
pcond_factor = NaN*ones(1,n);
pcond_factor_err = NaN*ones(1,n);
time_deriv = NaN*ones(1,n);
time_deriv_err = NaN*ones(1,n);
sta_mean = NaN*ones(1,n);
sta_rms = NaN*ones(1,n);
sta_SAL = NaN*ones(m,n);
sta_SAL1 = NaN*ones(m,n);
sta_SAL_err = NaN*ones(m,n);
sta_COND = NaN*ones(m,n);
sta_COND_err = NaN*ones(m,n);
fcoef = [];
fbreaks = [];


% -------------------------------------
% STEP 1 - Theta calseries segments
% -------------------------------------

% compute the number of segments for theta calseries
unique_cal = unique(calseries);
bad = find(unique_cal == 0); % 0 denotes bad profile to be skipped
if ~isempty(bad)
    unique_cal(bad) = [];
end
n_seq = length(unique_cal);
fprintf(' INFO: %d segment(s) used for the selection of theta levels\n',n_seq)

% Initialisation
Theta=NaN*zeros(10,n_seq);
Index=NaN*zeros(10,n);
Plevels=NaN*zeros(10,n_seq);

% loop through sequences of theta calseries
for iseq=1:n_seq

    calindex = find(calseries==unique_cal(iseq));
    k = length(calindex);

    unique_SAL = SAL(:, calindex);
    unique_PTMP = PTMP(:, calindex);
    unique_PRES = PRES(:, calindex);
    unique_la_ptmp = la_ptmp(:, calindex);

    % DD (2026/06/30) correct Theta by Theta(:,iseq)
    [Theta(:,iseq), P, index, var_s_th, th] =...
     find_10thetas( unique_SAL, unique_PTMP, unique_PRES, unique_la_ptmp, use_theta_gt, use_theta_lt, use_pres_gt, use_pres_lt, use_percent_gt);

    % DD (2024/09-4): Save Theta-related information for plotting functions
    Index(:,calindex)=index;
    Plevels(:,iseq)=P;
    [mvth,~]=size(var_s_th);
    if iseq==1
        Var_s_Thetas=NaN.*ones(mvth,n_seq);
        Thetas=NaN.*ones(mvth,n_seq);
        Var_s_Thetas(:,iseq)=var_s_th;
        Thetas(:,iseq)=th;
    else
        [mVth,~]=size(Var_s_Thetas);
        tmp1=Var_s_Thetas;
        tmp2=Thetas;
        mm=max(mVth,mvth);
        Var_s_Thetas=NaN.*ones(mm,n_seq);
        Thetas=NaN.*ones(mm,n_seq);
        Var_s_Thetas(1:mVth,:)=tmp1;
        Var_s_Thetas(1:mvth,iseq)=var_s_th;
        Thetas(1:mVth,:)=tmp2;
        Thetas(1:mvth,iseq)=th;
    end


end

ls_theta = strcat( po_system_configuration.FLOAT_CALIB_DIRECTORY, ...
        pn_float_dir, 'selected_theta_', pn_float_name,'.mat');
save(ls_theta, 'Theta','Index','Plevels','Var_s_Thetas','Thetas')

% ---------------------------------------------------
% STEP 1 - Continuous piecewisefit calseries segments
% ---------------------------------------------------

% compute the number of segments for continuous piecewisefit calseries
sstatus = 1;
unique_cal = unique(calseries_for_cont_pwf);
bad = find(unique_cal == 0); % 0 denotes bad profile to be skipped
if ~isempty(bad)
    unique_cal(bad) = [];
end
n_seq = length(unique_cal);
fprintf(' INFO: %d segment(s) used for continuous piecewisefit calibration\n',n_seq)

%Check max_break configuration
l_max_breaks=length(max_breaks);
if n_seq == 1 & l_max_breaks > 1
    fprintf(' ERROR in specifying the maximum number of break points within max_breaks:\n the max_breaks array has dimension %d \n   it should be 1\n', ...
        l_max_breaks);
    sstatus = 0;
elseif n_seq > 1 % we have multiple cal series, make sure that break information is provided for all segments
    if l_max_breaks == 1 % only one max_break specified, specify max_breaks for all cal series segments
        max_breaks = ones(n_seq,1)*max_breaks;
    elseif l_max_breaks ~= n_seq % error in specification of max_breaks
        fprintf(' ERROR in specifying the maximum number of break points within max_breaks: \n   the max_breaks array has dimension %d \n   it should be either 1 or the number of piecewisefit segments: %d \n', ...
            [l_max_breaks n_seq]);
        sstatus = 0;
    end
end

if sstatus == 0 % set_calseries returned a bad status variable, go ahead and write out file with NaNs
    ls_float_calib_filename = strcat( po_system_configuration.FLOAT_CALIB_DIRECTORY, pn_float_dir, po_system_configuration.FLOAT_CALIB_PREFIX, pn_float_name, po_system_configuration.FLOAT_CALIB_POSTFIX ) ;

    save( ls_float_calib_filename, 'cal_SAL', 'cal_SAL_err', 'pcond_factor', ...
        'pcond_factor_err', 'cal_COND', 'cal_COND_err', 'time_deriv', 'time_deriv_err', ...
          'sta_mean', 'sta_rms', 'sta_SAL', 'sta_SAL_err', 'PROFILE_NO', 'fcoef', 'fbreaks' )
    return
end


% Check breaks configuration
if ~isempty(breaks)
    breaks = breaks((isfinite(breaks))); % security
    [n1, n2] = size(breaks);
    nb=max(n1,n2);
    fprintf(' INFO: %d fixed break(s) prescribed\n', nb)
    
    if ~isempty(max_breaks)
        if n_seq ==1
            if nb > max_breaks
                max_breaks = nb;
                fprintf(' WARNING: max_breaks value was lower than the number of breaks indicated in breaks, max_breaks was overwritten with %d\n',nb)
            end
        elseif n_seq>1
            for i_seq=1:n_seq
                calindex = find(calseries_for_cont_pwf==unique_cal(i_seq));
                breaks_iseq=intersect(breaks, calindex);
                if ~isempty(breaks_iseq)
                    [n1, n2] = size(breaks_iseq);
                    nb_iseq=max(n1,n2);
                    fprintf('  INFO %d fixed break(s) prescribed for segment %d\n', nb_iseq,i_seq)
                    if nb_iseq>max_breaks(i_seq)
                        max_breaks(i_seq)=nb_iseq;
                        fprintf('   WARNING: max_breaks value for segment %d was lower than the number of breaks indicated in breaks, max_breaks was overwritten with %d\n',i_seq,nb_iseq)
                    end
                else
                    fprintf('  INFO no fixed break(s) prescribed for segment %d\n', i_seq)
                end
            end
        end
    end
end



% loop through sequences of continuous piecewisefit calseries segments
for iseq=1:n_seq
    
    fprintf(' INFO: dealing with segment %d for continuous piecewisefit calibration \n',iseq)

    calindex = find(calseries_for_cont_pwf==unique_cal(iseq));
    k = length(calindex);
    index=Index(:,calindex);
    unique_SAL = SAL(:,calindex);
    unique_PTMP = PTMP(:,calindex);
    unique_PRES = PRES(:,calindex);
    unique_mapped_sal = mapped_sal(:,calindex);
    unique_mapsalerrors = mapsalerrors(:,calindex);
    unique_coord_float = coord_float(calindex,:); 
    
    %initialisation
    ten_SAL = NaN.*ones(10,k);
    ten_PTMP = NaN.*ones(10,k);
    ten_PRES = NaN.*ones(10,k);
    ten_mapped_sal = NaN.*ones(10,k);
    ten_mapsalerrors = NaN.*ones(10,k);
    

    pp=find(isnan(index)==0);

    if(isempty(pp)==0) % only proceed when there are valid levels ----
        
        for ipr=1:k
            jj=find(isnan(index(:,ipr))==0);
            if(isempty(jj)==0)
                ith=index(jj,ipr);
                ten_SAL(1:length(jj),ipr) = unique_SAL(ith, ipr);
                ten_PTMP(1:length(jj),ipr) = unique_PTMP(ith, ipr);
                ten_PRES(1:length(jj),ipr) = unique_PRES(ith, ipr);
                ten_mapped_sal(1:length(jj),ipr) = unique_mapped_sal(ith, ipr);
                ten_mapsalerrors(1:length(jj),ipr) = unique_mapsalerrors(ith, ipr);
            end
        end

        % calculate potential conductivities and errors for mapped values and float values
        % calculate pcond error by perturbing salinity ... to avoid problems caused by non-linearity of the Equation of State ---

        ICOND = sw_c3515*sw_cndr(ten_SAL, ten_PTMP, 0);
        mapped_cond = sw_c3515*sw_cndr(ten_mapped_sal, ten_PTMP, 0);
        mapped_cond1 = sw_c3515*sw_cndr(ten_mapped_sal+ten_mapsalerrors/100, ten_PTMP, 0);
        mapconderrors = 100*abs(mapped_cond-mapped_cond1);

        x = x_in(:,calindex); % independent variable for piecewise fit (Profile Number)
        y = mapped_cond./ICOND; % dependent variable for fit (conductivity ratio)
        err = mapconderrors./ICOND; % error estimate for dependent variable (in ratio form)

        % calculate off-diagonal terms for error estimate --------
        %covariance = build_ptmp_cov(ten_PTMP); % build the data covariance matrix   % vertical covariance only
        %covariance = build_ptmp_xy_cov(ten_PTMP,unique_coord_float,po_system_configuration);   % change config 129 :  vertical and horizontal covariances
        covariance = build_ptmp_xyt_cov(ten_PTMP,unique_coord_float,po_system_configuration);   % ccabanes 01/06/2020


        % for debugging purposes to speed up calculations, use next line for first time calculation
        % and then comment out the call to build_ptmp_cov and load the covariance matrix
        %###    eval(['save ' strcat( po_system_configuration.FLOAT_CALIB_DIRECTORY, pn_float_dir, po_system_configuration.FLOAT_CALSERIES_PREFIX , pn_float_name, 'cov.mat') ' covariance']);
        %###    eval(['load ' strcat( po_system_configuration.FLOAT_CALIB_DIRECTORY, pn_float_dir, po_system_configuration.FLOAT_CALSERIES_PREFIX , pn_float_name, 'cov.mat') ' covariance']);
        % for debugging

        % use covariance to estimate off diagonal error terms
        % i.e. have weighting matrix include off diagonal terms -----

        % selecting the prescribed breaks for the current segment iseq
        breaks_iseq=intersect(breaks, calindex);
        
        % if no breaks points are set
        if isempty(breaks_iseq)
            [xfit(calindex), pcond_factor(calindex), pcond_factor_err(calindex), time_deriv(calindex), ...
             time_deriv_err(calindex), sta_mean(calindex), sta_rms(calindex), NDF, fitcoef, fitbreaks] = ...
              fit_cond(x, y, err, covariance, 'max_no_breaks', max_breaks(iseq));
        else
            if isempty(max_breaks(iseq))
                [xfit(calindex), pcond_factor(calindex), pcond_factor_err(calindex), time_deriv(calindex), ...
                 time_deriv_err(calindex), sta_mean(calindex), sta_rms(calindex), NDF, fitcoef, fitbreaks] = ...
                  fit_cond(x, y, err, covariance, 'breaks', breaks_iseq);
            else
                [xfit(calindex), pcond_factor(calindex), pcond_factor_err(calindex), time_deriv(calindex), ...
                 time_deriv_err(calindex), sta_mean(calindex), sta_rms(calindex), NDF, fitcoef, fitbreaks] = ...
                  fit_cond(x, y, err, covariance, 'breaks', breaks_iseq, 'max_no_breaks', max_breaks(iseq)); 
            end
        end

        % debug
        % ls_float_calib_dbg_filename = strcat( po_system_configuration.FLOAT_CALIB_DIRECTORY, pn_float_dir, ...
        %   po_system_configuration.FLOAT_CALIB_PREFIX, pn_float_name, "_" , num2str(iseq), "_b", po_system_configuration.FLOAT_CALIB_POSTFIX ) ;
        % save( ls_float_calib_dbg_filename,'x','y','err','covariance','calindex','ten_mapped_sal','ten_PTMP','ten_SAL')
        % end debug


        % apply calibrations to float data ------------------------

        if ~isempty(pcond_factor(calindex))
            unique_COND = sw_c3515*sw_cndr( unique_SAL, unique_PTMP, 0);
            cal_COND(:,calindex) = ( ones(m,1)*pcond_factor(calindex) ).*unique_COND;
            cal_SAL(:,calindex) = sw_salt( cal_COND(:,calindex)/sw_c3515, unique_PTMP, 0);

            % estimate the error in salinity ---------------------------------
            cal_COND_err(:,calindex) = ( ones(m,1)*pcond_factor_err(calindex) ).*unique_COND;
            cal_SAL1(:,calindex) = sw_salt( (cal_COND(:,calindex)+cal_COND_err(:,calindex))/sw_c3515, unique_PTMP, 0 );
            cal_SAL_err(:,calindex) = abs(cal_SAL(:,calindex)-cal_SAL1(:,calindex));

            % estimate the error in salinity for station by station fit ----
            sta_COND(:,calindex) = ( ones(m,1)*sta_mean(calindex) ).*unique_COND;
            sta_SAL(:,calindex) = sw_salt( sta_COND(:,calindex)/sw_c3515, unique_PTMP, 0);
            sta_COND_err(:,calindex) = ( ones(m,1)*sta_rms(calindex) ).*unique_COND;
            sta_SAL1(:,calindex) = sw_salt( (sta_COND(:,calindex)+sta_COND_err(:,calindex))/sw_c3515, unique_PTMP, 0 );
            sta_SAL_err(:,calindex) = abs(sta_SAL(:,calindex)-sta_SAL1(:,calindex));

            fcoef(iseq,1:length(fitcoef)) = fitcoef;
            if ~isempty(fitbreaks)
                fbreaks(iseq,1:length(fitbreaks)) = fitbreaks;
            end

        end

    end %if there are valid levels
end %for each unique_cal



% save calibration data --------------------------------

ls_float_calib_filename = strcat( po_system_configuration.FLOAT_CALIB_DIRECTORY, pn_float_dir, po_system_configuration.FLOAT_CALIB_PREFIX, pn_float_name, po_system_configuration.FLOAT_CALIB_POSTFIX ) ;

%save( ls_float_calib_filename, 'cal_SAL', 'cal_SAL_err', 'pcond_factor', 'pcond_factor_err', 'cal_COND', 'cal_COND_err', ...
%       'time_deriv', 'time_deriv_err', 'sta_mean', 'sta_rms', 'sta_SAL', 'sta_SAL_err', 'PROFILE_NO', 'fcoef', 'fbreaks' )

save( ls_float_calib_filename, 'cal_SAL', 'cal_SAL_err', 'pcond_factor', 'pcond_factor_err',  ...
       'time_deriv', 'time_deriv_err', 'sta_mean', 'sta_rms', 'sta_SAL', 'sta_SAL_err', 'PROFILE_NO', 'fcoef', 'fbreaks' )
