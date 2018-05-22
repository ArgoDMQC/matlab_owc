 
function calculate_piecewisefit( pn_float_dir, pn_float_name, po_system_configuration )

%
% Annie Wong, October 2008
% Breck Owens, November 2007
%

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


% load calibration settings -----------------

lo_float_calseries = load( strcat( po_system_configuration.FLOAT_CALIB_DIRECTORY, pn_float_dir, po_system_configuration.FLOAT_CALSERIES_PREFIX , pn_float_name, po_system_configuration.FLOAT_MAPPED_POSTFIX  ) );

calseries = lo_float_calseries.calseries;
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
cal_SAL_err = NaN*ones(m,n);
cal_COND = NaN*ones(m,n);
cal_COND_err = NaN*ones(m,n);
pcond_factor = NaN*ones(1,n);
pcond_factor_err = NaN*ones(1,n);
time_deriv = NaN*ones(1,n);
time_deriv_err = NaN*ones(1,n);
sta_mean = NaN*ones(1,n);
sta_rms = NaN*ones(1,n);
sta_SAL = NaN*ones(m,n);
sta_SAL_err = NaN*ones(m,n);
fcoef = [];
fbreaks = [];


% for each unique calseries -----------

sstatus = 1;
unique_cal = unique(calseries);
bad = find(unique_cal == 0); % 0 denotes bad profile to be skipped
if ~isempty(bad)
    unique_cal(bad) = [];
end
n_seq = length(unique_cal);

if n_seq == 1 & length(max_breaks) > 1
    disp(sprintf(' ERROR in specifying the number of possible break points\n %d max_breaks specified, should be %d value', ...
        [length(max_breaks) n_seq]) );
    sstatus = 0;
elseif n_seq > 1 % we have multiple cal series, make sure that break information is provided for all segments
    if length(max_breaks) == 1 % only one max_break specified, specify max_breaks for all cal series segments
        max_breaks = ones(n_seq,1)*max_breaks;
    elseif length(max_breaks) ~= n_seq % error in specification of max_breaks
        disp(sprintf(' ERROR in specifying the number of possible break points\n %d max_breaks specified, should be either 1 or %d values', ...
            [length(max_breaks) n_seq]) );
        sstatus = 0;
    end
end

if ~isempty(breaks)
    [ns nb] = size(breaks);
    if ns ~= n_seq % error in specifying breaks
        disp(sprintf(' ERROR in specifying break points\n For multiple cal series, need to specify breaks for each series\n Have %d cal series and %d sets of breaks', ...
            [n_seq ns]) );
        sstatus = 0;
    end
    for n=1:n_seq
        nb = length(find(isfinite(breaks(n,:))));
        if nb > max_breaks(n)
            disp(sprintf(' ERROR, for cal series %d max number of breaks %d less than %d prescribed breaks',unique_cal(n),max_breaks(n),nb))
            sstatus = 0;
        elseif nb < max_breaks(n)
            disp(sprintf(' Specified %d breaks.  Will search for up to %d breaks.',nb,max_breaks(n)))
        else
            disp(sprintf('  %d fixed breaks prescribed', nb))
        end
    end
end

if sstatus == 0 % set_calseries returned a bad status variable, go ahead and write out file with NaNs
    ls_float_calib_filename = strcat( po_system_configuration.FLOAT_CALIB_DIRECTORY, pn_float_dir, po_system_configuration.FLOAT_CALIB_PREFIX, pn_float_name, po_system_configuration.FLOAT_CALIB_POSTFIX ) ;

    save( ls_float_calib_filename, 'cal_SAL', 'cal_SAL_err', 'pcond_factor', ...
        'pcond_factor_err', 'cal_COND', 'cal_COND_err', 'time_deriv', 'time_deriv_err', ...
          'sta_mean', 'sta_rms', 'sta_SAL', 'sta_SAL_err', 'PROFILE_NO', 'fcoef', 'fbreaks' )
    return
end

% loop through sequences of calseries, ie if time series is broken into segments --
for i=1:n_seq

    calindex = find(calseries==unique_cal(i));
    k = length(calindex);

        % choose 10 float theta levels to use in the piecewise linear fit --------

        unique_SAL = SAL(:, calindex);
        unique_PTMP = PTMP(:, calindex);
        unique_PRES = PRES(:, calindex);
        unique_la_ptmp = la_ptmp(:, calindex);
        unique_mapped_sal = mapped_sal(:, calindex);
        unique_mapsalerrors = mapsalerrors(:, calindex);

        ten_SAL = NaN.*ones(10,k);
        ten_PTMP = NaN.*ones(10,k);
        ten_PRES = NaN.*ones(10,k);
        ten_mapped_sal = NaN.*ones(10,k);
        ten_mapsalerrors = NaN.*ones(10,k);

        [Theta, P, index, var_s_th, th] =...
         find_10thetas( unique_SAL, unique_PTMP, unique_PRES, unique_la_ptmp, use_theta_gt, use_theta_lt, use_pres_gt, use_pres_lt, use_percent_gt);

        pp=find(isnan(index)==0);
        if(isempty(pp)==0) % only proceed when there are valid levels ----

            for ipr=1:k
                jj=find(isnan(index(:,ipr))==0);
                if(isempty(jj)==0)
                    ten_SAL(1:length(jj),ipr) = unique_SAL( index(jj,ipr), ipr);
                    ten_PTMP(1:length(jj),ipr) = unique_PTMP( index(jj,ipr), ipr);
                    ten_PRES(1:length(jj),ipr) = unique_PRES( index(jj,ipr), ipr);
                    ten_mapped_sal(1:length(jj),ipr) = unique_mapped_sal( index(jj,ipr), ipr);
                    ten_mapsalerrors(1:length(jj),ipr) = unique_mapsalerrors( index(jj,ipr), ipr);
                end
            end

            % calculate potential conductivities and errors for mapped values and float values
            % calculate pcond error by perturbing salinity ... to avoid problems caused by non-linearity of the Equation of State ---

            ICOND = sw_c3515*sw_cndr( ten_SAL, ten_PTMP, 0);
            mapped_cond = sw_c3515*sw_cndr( ten_mapped_sal, ten_PTMP, 0);
            mapped_cond1 = sw_c3515*sw_cndr( ten_mapped_sal+ten_mapsalerrors/100, ten_PTMP, 0);
            mapconderrors = 100*abs(mapped_cond-mapped_cond1);

            x = x_in(:,calindex); % independent variable for piecewise fit (Profile Number)
            y = mapped_cond./ICOND; % dependent variable for fit (conductivity ratio)
            err = mapconderrors./ICOND; % error estimate for dependent variable (in ratio form)

            % calculate off-diagonal terms for error estimate --------

            covariance = build_ptmp_cov(ten_PTMP); % build the data covariance matrix

            % for debugging purposes to speed up calculations, use next line for first time calculation
            % and then comment out the call to build_ptmp_cov and load the covariance matrix

            %###    eval(['save ' strcat( po_system_configuration.FLOAT_CALIB_DIRECTORY, pn_float_dir, po_system_configuration.FLOAT_CALSERIES_PREFIX , pn_float_name, 'cov.mat') ' covariance']);
            %###    eval(['load ' strcat( po_system_configuration.FLOAT_CALIB_DIRECTORY, pn_float_dir, po_system_configuration.FLOAT_CALSERIES_PREFIX , pn_float_name, 'cov.mat') ' covariance']);
            % for debugging

            % use covariance to estimate off diagonal error terms
            % i.e. have weighting matrix include off diagonal terms -----

            % if no breaks points are set
            if isempty(breaks)
                [xfit(calindex), pcond_factor(calindex), pcond_factor_err(calindex), time_deriv(calindex), ...
                 time_deriv_err(calindex), sta_mean(calindex), sta_rms(calindex), NDF, fitcoef, fitbreaks] = ...
                  fit_cond(x, y, err, covariance, 'max_no_breaks', max_breaks(i));
            else
                breaks_in = breaks(i,:);
                breaks_in = breaks_in(find(isfinite(breaks_in)));
                if isempty(max_breaks(i))
                    [xfit(calindex), pcond_factor(calindex), pcond_factor_err(calindex), time_deriv(calindex), ...
                     time_deriv_err(calindex), sta_mean(calindex), sta_rms(calindex), NDF, fitcoef, fitbreaks] = ...
                      fit_cond(x, y, err, covariance, 'breaks', breaks_in);
                else
                    [xfit(calindex), pcond_factor(calindex), pcond_factor_err(calindex), time_deriv(calindex), ...
                     time_deriv_err(calindex), sta_mean(calindex), sta_rms(calindex), NDF, fitcoef, fitbreaks] = ...
                      fit_cond(x, y, err, covariance, 'breaks', breaks_in, 'max_no_breaks', max_breaks(i));
                end
            end


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

                fcoef(i,1:length(fitcoef)) = fitcoef;
                if ~isempty(fitbreaks)
                    fbreaks(i,1:length(fitbreaks)) = fitbreaks;
                end

            end

        end %if there are valid levels
end %for each unique_cal


% save calibration data --------------------------------

ls_float_calib_filename = strcat( po_system_configuration.FLOAT_CALIB_DIRECTORY, pn_float_dir, po_system_configuration.FLOAT_CALIB_PREFIX, pn_float_name, po_system_configuration.FLOAT_CALIB_POSTFIX ) ;

save( ls_float_calib_filename, 'cal_SAL', 'cal_SAL_err', 'pcond_factor', 'pcond_factor_err', 'cal_COND', 'cal_COND_err', ...
       'time_deriv', 'time_deriv_err', 'sta_mean', 'sta_rms', 'sta_SAL', 'sta_SAL_err', 'PROFILE_NO', 'fcoef', 'fbreaks' )

