% Created:
%   [DD/MM/YYYY]
%   [??/10/2006] Breck Owens
%
% Updated:
%   [21/05/2025] Dirk Slawinski
%   - [DS2025] add flags to toggle headless mode
%   [??/??/2023] Cecile Cabanes / cc
%   - add split
%   [26/07/2021] Dirk Slawinski
%   - set plotting to be 'headless' and append 'argo' and 'ctd' depending
%   on climatology DB used
%   [04/02/2021] JLL
%   - change output to PNG not EPS
%   [4/6/2011] Annie Wong
%   -  adjusted
%   [??/10/2006] Breck Owens
%   - created file

% [DS2025/]
% - add goHeadless flag
function plot_diagnostics_ow_CSIRO( pn_float_dir, pn_float_name, po_system_configuration , ...
    goHeadless)
% [/DS2025]

%   pn_float_dir='uw/';
%   pn_float_name='R5902134';
%   po_system_configuration = load_configuration( 'ow_config.txt' );

    close all
    
% [DS2025/]
    if(~exist("goHeadless", "var"))
        goHeadless = false;
    end % if(~exist("goHeadless"))

    % set the plot type, -dpng vs -depsc
    pltFileType = '-dpng';
% [/DS2025]



    % modify pn_float_name for title if the name contains '_' ------------
    ii = find(pn_float_name == '_');
    if(isempty(ii) == 0)
        title_floatname = strcat( pn_float_name(1:ii - 1), ...
            '\_', pn_float_name(ii + 1:length(pn_float_name)));
    else % if(isempty(ii) == 0)
        title_floatname = pn_float_name;
    end % if(isempty(ii) == 0)


    % load data from /float_source, /float_mapped, /float_calib, and config --------------

    lo_float_source_data = load(fullfile(po_system_configuration.FLOAT_SOURCE_DIRECTORY, ...
        pn_float_dir, strcat( pn_float_name, po_system_configuration.FLOAT_SOURCE_POSTFIX ) ) );

    PROFILE_NO = lo_float_source_data.PROFILE_NO;
    LAT  = lo_float_source_data.LAT;
    LONG = lo_float_source_data.LONG;
    PRES = lo_float_source_data.PRES;
    TEMP = lo_float_source_data.TEMP;
    PTMP = lo_float_source_data.PTMP;
    SAL  = lo_float_source_data.SAL;

    if(isempty(find(isnan(PRES) == 0)) == 0)
        % if no data exists, terminate here, no plots will be produced

        lo_float_mapped_data = load( fullfile(po_system_configuration.FLOAT_MAPPED_DIRECTORY, ...
            pn_float_dir, strcat( po_system_configuration.FLOAT_MAPPED_PREFIX, pn_float_name, ...
            po_system_configuration.FLOAT_MAPPED_POSTFIX ) ) ) ;

        mapped_sal  = lo_float_mapped_data.la_mapped_sal;
        mapsalerrors = lo_float_mapped_data.la_mapsalerrors;
        la_ptmp = lo_float_mapped_data.la_ptmp;
        selected_hist = lo_float_mapped_data.selected_hist;

        if(isempty(find(isnan(mapped_sal) == 0)) == 0)
            % if mapping exists, terminate here, no plots will be produced

            lo_float_calib_data = load(fullfile(po_system_configuration.FLOAT_CALIB_DIRECTORY, ...
                pn_float_dir, strcat( po_system_configuration.FLOAT_CALIB_PREFIX, pn_float_name, ...
                po_system_configuration.FLOAT_CALIB_POSTFIX ) ) );

            cal_SAL = lo_float_calib_data.cal_SAL;
            cal_SAL_err = lo_float_calib_data.cal_SAL_err;
            pcond_factor = lo_float_calib_data.pcond_factor;
            pcond_factor_err = lo_float_calib_data.pcond_factor_err;

            lo_float_calseries = load( fullfile(po_system_configuration.FLOAT_CALIB_DIRECTORY, ...
                pn_float_dir, strcat( po_system_configuration.FLOAT_CALSERIES_PREFIX, pn_float_name, ...
                po_system_configuration.FLOAT_MAPPED_POSTFIX ) ) );

            use_theta_gt = lo_float_calseries.use_theta_gt;
            use_theta_lt = lo_float_calseries.use_theta_lt;
            use_pres_gt = lo_float_calseries.use_pres_gt;
            use_pres_lt = lo_float_calseries.use_pres_lt;
            use_percent_gt = lo_float_calseries.use_percent_gt;

            % load the station by station fits
            load(fullfile( ...
                po_system_configuration.FLOAT_CALIB_DIRECTORY, pn_float_dir, strcat( ...
                po_system_configuration.FLOAT_CALIB_PREFIX, pn_float_name, ...
                po_system_configuration.FLOAT_CALIB_POSTFIX ) ), ...
                '-regexp', '^sta');


            % plot the float locations (figure 1) -----------------------

            load( fullfile( po_system_configuration.CONFIG_DIRECTORY, ...
                po_system_configuration.CONFIG_COASTLINES ), ...
                'coastdata_x', 'coastdata_y' );

            [m,n] = size(PRES);

% [DS2025/]
%            figure;
            if(goHeadless)
                figure('Visible', 'off');
            else % if(goHeadless)
                figure('Visible', 'on');
            end % if(goHeadless)
% [/DS2025]
            set(gcf, 'defaultaxeslinewidth', 2)
            set(gcf, 'defaultlinelinewidth', 2)
            set(gcf, 'defaultaxesfontsize', 16)
            
            colormap(jet(n));
            c = colormap;

            x = [];
            y = [];
            if(isempty(selected_hist) == 0)
              x = selected_hist(:, 1);
              y = selected_hist(:, 2);
            end % if(isempty(selected_hist) == 0)

            if( max(LONG)+30 > 360 | min(LONG)-30 < 0 ) 
                % if there are float data within 30 deg lon of the prime meridian
                if( isempty( find(LAT < -50 & LONG > 250 & LONG < 285) ) == 0)
                    jj = find(LONG < LONG(1) - 30); % and if there are float data west of the Drake Passage
                    LONG(jj) = LONG(jj) + 360;    % then use the first float location point minus 30 deg lon
                    kk = find(x < LONG(1) - 30);    % as the western edge of the map
                    x(kk) = x(kk) + 360;
                    mm = find(coastdata_x < LONG(1) - 30);
                    coastdata_x(mm) = coastdata_x(mm) + 360;
                    nn = find(coastdata_x < LONG(1) - 30 + 1 | coastdata_x > LONG(1) - 30 + 360 - 1); %remove coastline wraparound
                    coastdata_x(nn) = NaN;
                else % if( isempty( find(LAT<-50 & LONG>250 & LONG<285) ) == 0)
                    jj = find(LONG < 250);     % or else if there are no float data west of the Drake Passage
                    LONG(jj) = LONG(jj) + 360; % then use 250 deg lon as the western edge of the map
                    kk = find(x < 250);
                    x(kk) = x(kk) + 360;
                    mm = find(coastdata_x < 250);
                    coastdata_x(mm) = coastdata_x(mm) + 360;
                    nn = find(coastdata_x < 251 | coastdata_x > 609); %remove coastline wraparound
                    coastdata_x(nn) = NaN;
                end % if( isempty( find(LAT<-50 & LONG>250 & LONG<285) ) == 0)
            end % if( max(LONG)+30 > 360 | min(LONG)-30 < 0 ) 

            % float path
            plot(LONG, LAT, 'r-');
            hold on
            % historical points
            plot(x, y, 'b.')
% [DS2025/]
%            legend('float','historical points', 'Location', 'Best')
%            % coast line data
%            plot(coastdata_x, coastdata_y, 'k.-');

%            for i = 1:n
%                h = plot(LONG(i), LAT(i), '+');
%                set(h, 'color', c(i, :));
%                j = text(LONG(i), LAT(i), int2str(PROFILE_NO(i)));
%                set(j, 'color', c(i,:),'fontsize', 12, 'hor', 'cen');
%            end % for i = 1:n
            LLscat = scatter(LONG, LAT, '+');
            LLscat.CData = c;
            LLtxt = text(LONG, LAT, strsplit(num2str(PROFILE_NO(:)')));
            for i = 1:n
                LLtxt(i).Color = c(i, :);
            end % for i = 1:n
            % coast line data
            plot(coastdata_x, coastdata_y, 'k.-');
            legend('float', 'historical points', 'location','Best')
% [/DS2025]
            axis([min(LONG) - 30, max(LONG) + 30, ...
                min(LAT) - 20, max(LAT) + 20]);
            set(gca, 'FontSize', 12);
            xlabel('Longitude');
            ylabel('Latitude');
% [DS2025/]
%            title( strcat(title_floatname, ' profile locations with historical data' ) );
            title([title_floatname, ' profile locations with historical data']);
% [/DS2025]

            h = gca;
            xticks = get(h, 'XTick');
            xticklabels = get(h, 'XTickLabel');
            ii = find(xticks > 360);
            xticks(ii) = xticks(ii) - 360;

            for j = 1:length(ii)
                enough = num2str(ones(3, 1));
                c = char(num2str(xticks(ii(j))), enough');
                if ischar(xticklabels)
                    xticklabels(ii(j), :) = c(1, :);
                end % if ischar(xticklabels)
                if iscell(xticklabels)
                    xticklabels(ii(j), :) = cellstr(c(1, :));  % Matt Donnelly: fix issue in Matlab 2016b
                end % if iscell(xticklabels)
            end % for j = 1:length(ii)
            set(gca, 'XTickLabel', xticklabels);

            drawnow
% [DS2025/]
% - matlab indicates that 'paperorientation' is no longer a valid parameter
%            set(gcf, 'papertype', 'usletter', 'paperunits', 'inches', 'paperorientation', ...
%                'portrait', 'paperposition', [0.25, 0.75, 8.0, 9.5]);
            set(gcf, 'papertype', 'usletter', 'paperunits', 'inches', ...
                'paperposition', [0.25, 0.75, 8.0, 9.5]);
% [/DS2025]

% [DS2025/]
%            print('-depsc', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, ...
%                pn_float_dir, pn_float_name, '_1.eps'));
            hasCtd = strfind(po_system_configuration.CONFIG_WMO_BOXES, ...
                'ctd');
            if(~isempty(hasCtd))
                fname = strcat(...
                    po_system_configuration.FLOAT_PLOTS_DIRECTORY, ...
                    pn_float_dir, pn_float_name, '_1_ctd');
            else % if(~isempty(hasCtd))
                fname = strcat(...
                    po_system_configuration.FLOAT_PLOTS_DIRECTORY, ...
                    pn_float_dir, pn_float_name, '_1');
            end % if(~isempty(hasCtd))
            print(pltFileType, '-r300', fname);
% [/DS2025]

            % plot the uncalibrated theta-S curves from the float (figure 2) --------
% [DS2025/]
%            figure;
            if(goHeadless)
                figure('Visible', 'off');
            else
                figure('Visible', 'on');
            end % if(goHeadless)
% [/DS2025]
            set(gcf, 'defaultaxeslinewidth', 2)
            set(gcf, 'defaultlinelinewidth', 2)
            set(gcf, 'defaultaxesfontsize', 14)
            
            jj = [1:ceil(n/30):n]; % the legend can only fit 30 profiles
            ok = [];
            for k = 1:length(jj)
                if ~isempty(find(isfinite(SAL(:, jj(k)))))
                    ok = [ok jj(k)];
                end % if ~isempty(find(isfinite(SAL(:, jj(k)))))
            end % for k = 1:length(jj)
            jj = ok;
            
            colormap(jet(n));
            c = colormap;
% [DS2025/]
%            qq = plot(PTMP(:,jj), SAL(:,jj));
            qq = plot(SAL(:, jj), PTMP(:, jj));
% [/DS2025]
            for i = 1:length(jj)
                set(qq(i), 'color', c(jj(i), :)) ;
            end % for i = 1:length(jj)
            hold on
% [DS2025/]
% - adding the legend should be the last thing to add if you are planning
% to limit what it shows as matla will just append what is added
%            legend(int2str([PROFILE_NO(jj)]'), ...
%                'Location', 'NorthEastOutside')
% [/DS2025]

            for i = 1:n
                % plot all remaining profiles
% [DS2025/] 
%               qq = plot( PTMP(:,i), SAL(:,i) );
                qq = plot(SAL(:, i), PTMP(:, i));
% [/DS2025]
                set(qq, 'color', c(i, :)) ;
            end % for i = 1:n

            [tlevels, plevels, index, var_s_Thetalevels, Thetalevels] = find_10thetas( ...
                SAL, PTMP, PRES, la_ptmp, use_theta_gt, use_theta_lt, use_pres_gt, ...
                use_pres_lt, use_percent_gt);

            for i = 1:n
                % DS NOTE: couldn't we just use b = find(~isnan(... ?
                b = find( isnan(index(:, i)) == 0 );
                a = index(b, i);
                if ~isempty(a)
% [DS2025/]
%                    h = errorbar( la_ptmp(a, i), mapped_sal(a, i), mapsalerrors(a, i), ...
%                        'o', 'color', c(i, :) );
                    h = errorbar( mapped_sal(a, i), la_ptmp(a, i), mapsalerrors(a, i), 'horizontal', ...
                        'o', 'color', c(i, :));
% [/DS2025]
                end % if ~isempty(a)
            end % for i = 1:n
% [DS2025/]
            legend(int2str([PROFILE_NO(jj)]'), ...
                'Location', 'NorthEastOutside')
% [/DS2025]

% [DS2025/]
% - no need to rotate and create a 3d plot since errorbars can be horizontal
%            view([90 -90]);
% [/DS2025]
% 
            set(gca, 'FontSize', 12)
% [DS2025/]
%            ylabel('Salinity (PSS-78)');
%            xlabel('\theta ^{\circ} C');
            xlabel('Salinity (PSS-78)');
            ylabel('\theta ^{\circ} C');
% [/DS2025]

            max_s = max([max(SAL), max(mapped_sal)]) + 0.1;
            min_s = min([min(SAL), min(mapped_sal)]) - 0.1;
            max_t = max(max(PTMP)) + 1;
            min_t = min(min(PTMP)) - 1;
% [DS2025/]
%            axis([min_t, max_t, min_s, max_s]);
            axis([min_s, max_s, min_t, max_t]);
% [/DS2025]

            drawnow
% [DS2025/]
%            title( strcat( title_floatname, ...
%                ' uncalibrated float data (-) and mapped salinity (o) with objective errors' ) );
            title({[title_floatname ' uncalibrated float data (-) and'], ...
                'mapped salinity (o) with objective errors'});
% [/DS2025]

% [DS2025/]
% - matlab indicates that 'paperorientation' is no longer a valid parameter
%            set(gcf, 'papertype', 'usletter', 'paperunits', 'inches', 'paperorientation', ...
%                'portrait', 'paperposition', [0.25, 0.75, 8.0, 9.5]);
            set(gcf, 'papertype', 'usletter', 'paperunits', 'inches', ...
                'paperposition', [0.25, 0.75, 8.0, 9.5]);
% [/DS2025]

% [DS2025/]
%            print('-depsc', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, pn_float_dir, ...
%                pn_float_name, '_2.eps'));
            hasCtd = strfind(po_system_configuration.CONFIG_WMO_BOXES, ...
                'ctd');
            if(~isempty(hasCtd))
                fname = strcat(...
                    po_system_configuration.FLOAT_PLOTS_DIRECTORY, ...
                    pn_float_dir, pn_float_name, '_2_ctd');
            else % if(~isempty(hasCtd))
                fname = strcat(...
                    po_system_configuration.FLOAT_PLOTS_DIRECTORY, ...
                    pn_float_dir, pn_float_name, '_2');
            end % if(~isempty(hasCtd))
            print(pltFileType, '-r300', fname);
            %disp(fname)
% [/DS2025]


            % calibration curve (figure 3) --------------------------
            if(isempty(find(isnan(cal_SAL) == 0)) == 0) 
                % if no cal exists, terminate here, no plot will be produced
                Soffset = cal_SAL-SAL;
                avg_Soffset = NaN .* ones(1, n);
                avg_Soffset_err = NaN .* ones(1, n);
                Staoffset = sta_SAL - SAL;
                avg_Staoffset = NaN .* ones(1, n);
                avg_Staoffset_err = NaN .* ones(1, n);

                for i = 1:n
                    ii = [];
                    ii = find(isnan(Soffset(:, i)) == 0);
                    if ~isempty(ii)
                        avg_Soffset(i) = mean(Soffset(ii, i));
                        avg_Soffset_err(i) = mean(cal_SAL_err(ii, i));
                    else % if ~isempty(ii)
                        avg_Soffset(i) = NaN;
                        avg_Soffset_err(i) = NaN;
                    end % if ~isempty(ii)
                    ii = find(isnan(Staoffset(:, i)) == 0);
                    if ~isempty(ii)
                        avg_Staoffset(i) = mean(Staoffset(ii, i));
                        avg_Staoffset_err(i) = mean(sta_SAL_err(ii, i));
                    else % if ~isempty(ii)
                        avg_Staoffset(i) = NaN;
                        avg_Staoffset_err(i) = NaN;
                    end % if ~isempty(ii)
                end % for i = 1:n

% [DS2025/]
%               figure;
                if(goHeadless)
                    figure('Visible', 'off');
                else % if(goHeadless)
                    figure('Visible', 'on');
                end % if(goHeadless)
% [/DS2025]
                set(gcf, 'defaultaxeslinewidth', 2)
                set(gcf, 'defaultlinelinewidth', 2)
                set(gcf, 'defaultaxesfontsize', 16)
                subplot(2, 1, 1)
                plot(PROFILE_NO, pcond_factor, 'b-');
                hold on
                plot(PROFILE_NO, pcond_factor, 'g-');
                % plot station by station fit
                ok = find(isfinite(sta_mean));
                plot(PROFILE_NO(ok), sta_mean(ok), 'r-');
% [DS2025/]
%                legend('2 x cal error', '1 x cal error','1-1 profile fit', 'Location', 'Best');
% [/DS2205]
                errorbar(PROFILE_NO, pcond_factor, ...
                    2 * pcond_factor_err, 'b')
                errorbar(PROFILE_NO, pcond_factor, ...
                    pcond_factor_err, 'g*-')
                plot(PROFILE_NO(ok), sta_mean(ok), 'r-');

                plot( [0, max(PROFILE_NO)+1], [1,1], 'k-')
                axis([0, max(PROFILE_NO)+1, ...
                    min([pcond_factor-pcond_factor_err, 1]) - 0.0004, ...
                    max([pcond_factor+pcond_factor_err, 1]) + 0.0004 ])
                set(gca, 'FontSize', 12)
                ylabel('r') % multiplicative term has no units
% [DS2025/]
%                title( strcat(title_floatname, ...
%                    ' potential conductivity (mmho/cm) multiplicative correction r with errors'));
                title({[title_floatname 'potential conductivity (mmho/cm)'], ...
                    'multiplicative correction r with errors'})
                legend('2 x cal error', '1 x cal error','1-1 profile fit', 'Location', 'Best');
% [/DS2205]

                subplot(2, 1, 2)
                plot(PROFILE_NO, avg_Soffset, 'b-');
                hold on
                plot(PROFILE_NO, avg_Soffset, 'g-');
                % Plot station by station fit
                ok = find(isfinite(avg_Staoffset));
                plot(PROFILE_NO(ok), avg_Staoffset(ok), 'r-');
% [DS2025/]
%                legend('2 x cal error', '1 x cal error', '1-1 profile fit', 'Location', 'Best');
% [/DS2025]
                errorbar(PROFILE_NO, avg_Soffset, 2 * avg_Soffset_err, 'b')
                errorbar(PROFILE_NO, avg_Soffset, avg_Soffset_err, 'go-')
                plot(PROFILE_NO(ok), avg_Staoffset(ok), 'r-');
                
                axis([0, max(PROFILE_NO) + 1, ...
                    min([avg_Soffset - avg_Soffset_err, 0]) - 0.02, ...
                    max([avg_Soffset + avg_Soffset_err,0]) + 0.02 ])
                plot( [0, max(PROFILE_NO)+1], [0,0], 'k-')
                set(gca, 'FontSize', 12)
                xlabel('float profile number');
                ylabel('\Delta S')
% [DS2025/]
%                title( strcat(title_floatname, [' vertically-averaged salinity (PSS-78) ' ...
%                    'additive correction \Delta S with errors']) );
                title({[title_floatname ' vertically-averaged salinity (PSS-78) '], ...
                    'additive correction \Delta S with errors'})
                legend('2 x cal error', '1 x cal error', '1-1 profile fit', 'Location', 'Best');
% [/DS2025]

                drawnow
% [DS2025/]
% - matlab indicates that 'paperorientation' is no longer a valid parameter
%                set(gcf, 'papertype', 'usletter', 'paperunits', 'inches', 'paperorientation', ...
%                    'portrait',  'paperposition', [0.25, 0.75, 8.0, 9.5]);
                set(gcf, 'papertype', 'usletter', 'paperunits', 'inches', ...
                    'paperposition', [0.25, 0.75, 8.0, 9.5]);
% [/DS2025]

% [DS2025/]
%                 print('-depsc', strcat( po_system_configuration.FLOAT_PLOTS_DIRECTORY, ...
%                    pn_float_dir, pn_float_name, '_3.eps'));
                hasCtd = strfind(po_system_configuration.CONFIG_WMO_BOXES, ...
                    'ctd');
                if(~isempty(hasCtd))
                    fname = strcat(...
                        po_system_configuration.FLOAT_PLOTS_DIRECTORY, ...
                        pn_float_dir, pn_float_name, '_3_ctd');
                else % if(~isempty(hasCtd))
                    fname = strcat(...
                        po_system_configuration.FLOAT_PLOTS_DIRECTORY, ...
                        pn_float_dir, pn_float_name, '_3');
                end % if(~isempty(hasCtd))
                print(pltFileType, '-r300', fname);
% [/DS2025]

                % plot the calibrated theta-S curves from the float (figure 4) --------------------------
% [DS2025/]
%                figure;
                if(goHeadless)
                    figure('Visible', 'off');
                else % if(goHeadless)
                    figure('Visible', 'on');
                end % if(goHeadless)
% [/DS2025]
                set(gcf, 'defaultaxeslinewidth', 2)
                set(gcf, 'defaultlinelinewidth', 2)
                set(gcf, 'defaultaxesfontsize', 14)
                
                jj = [1:ceil(n/30):n]; % the legend can only fit 30 profiles
                ok = [];% check to make sure we have choosen profiles with good data
                for k = 1:length(jj)
                    if ~isempty(find(isfinite(SAL(:, jj(k)))))
                        ok = [ok jj(k)];
                    end % if ~isempty(find(isfinite(SAL(:, jj(k)))))
                end % for k = 1:length(jj)
                jj = ok;

                colormap(jet(n));
                c = colormap;
% [DS2025/]
%                qq = plot( PTMP(:, jj), cal_SAL(:, jj) );
                qq = plot(cal_SAL(:, jj), PTMP(:, jj));
% [/DS2025]
                for i = 1:length(jj)
                    set(qq(i), 'color', c(jj(i), :));
                end % for i = 1:length(jj)
                hold on;
% [DS2025/]
%                legend(int2str([PROFILE_NO(jj)]'), 'Location', 'NorthEastOutside')
% [/DS2025] 
               
                for i = 1:n
                    % plot all remaining profiles
% [DS2025/]
%                    qq = plot( PTMP(:, i), cal_SAL(:, i) );
                    qq = plot(cal_SAL(:, i), PTMP(:, i));
% [/DS2025]
                    set(qq, 'color', c(i, :)) ;
                end % for i = 1:n
                
                for i = 1:n
                    b = find( isnan(index(:,i))==0 );
                    a = index(b,i);
                    if ~isempty(a)
% [DS2025/]
%                        h = errorbar( la_ptmp(a, i), mapped_sal(a, i), mapsalerrors(a, i), ...
%                            'o', 'color', c(i, :));
                        h = errorbar(mapped_sal(a, i), la_ptmp(a, i), mapsalerrors(a, i), 'horizontal', ...
                            'o', 'color', c(i, :));
% [/DS2025]
                    end % if ~isempty(a)
                end % for i = 1:n

% [DS2025/]
%                view([90 -90]);
% [/DS2025]
                set(gca, 'FontSize', 12)

% [DS2025/]  
%                xlabel('\theta ^{\circ} C')
%                ylabel('Salinity (PSS-78)')
                ylabel('\theta ^{\circ} C')
                xlabel('Salinity (PSS-78)')
% [/DS2025]

                max_s = max([max(max(cal_SAL)), max(max(mapped_sal))]) + 0.1;
                min_s = min([min(min(cal_SAL)), min(min(mapped_sal))]) - 0.1;
                max_t = max(max(PTMP)) + 1;
                min_t = min(min(PTMP)) - 1;
                if(isnan(min_s) == 0)
% [DS2025/] 
%                    axis([min_t,max_t, min_s, max_s]);
                    axis([min_s,max_s, min_t, max_t]);
% [/DS2025]
                end % if(isnan(min_s) == 0)
                
                drawnow
% [DS2025/]
%                title( strcat(title_floatname, ...
%                    ' calibrated float data (-) and mapped salinity (o) with objective errors' ) );
                title({[title_floatname ' calibrated float data (-) and'], ...
                    'mapped salinity (o) with objective errors'});
                legend(int2str([PROFILE_NO(jj)]'), 'Location', 'NorthEastOutside')
% [/DS2025]
% [DS2025/]
% - matlab indicates that 'paperorientation' is no longer a valid parameter
%                set(gcf, 'papertype', 'usletter', 'paperunits', 'inches', 'paperorientation', ...
%                    'portrait', 'paperposition', [0.25, 0.75, 8.0 ,9.5]);
                set(gcf, 'papertype', 'usletter', 'paperunits', 'inches', ...
                    'paperposition', [0.25, 0.75, 8.0 ,9.5]);
% [/DS2025]

% [DS2025/]
%                print('-depsc', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, ...
%                    pn_float_dir, pn_float_name, '_4.eps'));
                hasCtd = strfind(po_system_configuration.CONFIG_WMO_BOXES, ...
                    'ctd');
                if(~isempty(hasCtd))
                    fname = strcat(...
                        po_system_configuration.FLOAT_PLOTS_DIRECTORY, ...
                        pn_float_dir, pn_float_name, '_4_ctd');
                else % if(~isempty(hasCtd))
                    fname = strcat(...
                        po_system_configuration.FLOAT_PLOTS_DIRECTORY, ...
                        pn_float_dir, pn_float_name, '_4');
                end % if(~isempty(hasCtd))
                print(pltFileType, '-r300', fname);
% [/DS2025]

                % Brian King's plot: salinity anomaly time series on theta levels (figure 5) ------------
% [DS2025/]
%                figure;
                if(goHeadless)
                    figure('Visible', 'off');
                else
                    figure('Visible', 'on');
                end % if(goHeadless)
% [/DS2025]
                set(gcf, 'defaultaxeslinewidth', 2)
                set(gcf, 'defaultlinelinewidth', 2)
                set(gcf, 'defaultaxesfontsize', 16)
                
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
                fl = anom(d, fl); % Brian King's routine
                subplot('position', [0.1 0.45 0.8 0.35])
% [DS2025/]
%                title(['       Salinity anom on theta.    ' title_floatname])
                title(['Salinity anom on theta.' title_floatname])
% [/DS2025]

                drawnow
% [DS2025/]
% - matlab indicates that 'paperorientation' is no longer a valid parameter
%                set(gcf, 'papertype', 'usletter', 'paperunits', 'inches', 'paperorientation', ...
%                    'portrait', 'paperposition', [0.25, 0.5, 8.0, 10.0]);
                set(gcf, 'papertype', 'usletter', 'paperunits', 'inches',  ...
                    'paperposition', [0.25, 0.5, 8.0, 10.0]);
% [/DS2025]

% [DS2025/]
%                print('-depsc', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, ...
%                    pn_float_dir, pn_float_name, '_5.eps'));
                hasCtd = strfind(po_system_configuration.CONFIG_WMO_BOXES, ...
                    'ctd');
                if(~isempty(hasCtd))
                    fname = strcat(...
                        po_system_configuration.FLOAT_PLOTS_DIRECTORY, ...
                        pn_float_dir, pn_float_name, '_5_ctd');
                else % if(~isempty(hasCtd))
                    fname = strcat(...
                        po_system_configuration.FLOAT_PLOTS_DIRECTORY, ...
                        pn_float_dir, pn_float_name, '_5');
                end % if(~isempty(hasCtd))
                print(pltFileType, '-r300', fname);
% [/DS2025]


                % plot salinity time series on theta levels with the smallest S variance (figure 6) ------------

                % CC changes 06/23 figure 6 to take into account the splitting of the time series 
                %  add l381 to l404
                unique_cal = unique(lo_float_calseries.calseries);
                n_seq = length(unique_cal);

                for iy = 1:n_seq
                    if unique_cal(iy) > 0
                        calindex = find(lo_float_calseries.calseries == unique_cal(iy));
                        k = length(calindex);
                        
                        % choose the two most stable theta level for each part of the time series --------
                        unique_SAL = SAL(:, calindex);
                        unique_PTMP = PTMP(:, calindex);
                        unique_PRES = PRES(:, calindex);
                        unique_la_ptmp = la_ptmp(:, calindex);
                        unique_mapped_sal = mapped_sal(:, calindex);
                        unique_mapsalerrors = mapsalerrors(:, calindex);
                        
                        ten_SAL = NaN .* ones(10, k);
                        ten_PTMP = NaN .* ones(10, k);
                        ten_PRES = NaN .* ones(10, k);
                        ten_mapped_sal = NaN .* ones(10, k);
                        ten_mapsalerrors = NaN .* ones(10, k);
                        
                        
                        [tlevels, plevels, index, var_s_Thetalevels, Thetalevels] = find_10thetas(...
                            unique_SAL, unique_PTMP, unique_PRES, unique_la_ptmp, ...
                            use_theta_gt, use_theta_lt, use_pres_gt, use_pres_lt, use_percent_gt);
                        %[tlevels, plevels, index, var_s_Thetalevels, Thetalevels] = find_10thetas( SAL, PTMP, PRES, la_ptmp, use_theta_gt, use_theta_lt, use_pres_gt, use_pres_lt, use_percent_gt); CC changes 06/23
                        
                        tplot = [1:2];
                        p = length(tplot);
                        Sint = NaN .* ones(p, n);
                        Smap = NaN .* ones(p, n);
                        Smaperr = NaN .* ones(p, n);
                        Scal = NaN .* ones(p, n);
                        Scalerr = NaN .* ones(p, n);
                        Thetalevel_indexes = NaN .* ones(p, n);
                        
                        % CC full time series is plotted on specified theta
                        trimPRES = PRES;  % use only manually specified THETA & PRES range ---
                        trimSAL = SAL;
                        trimPTMP = PTMP;
                        trim_mapped_sal = mapped_sal;
                        trim_mapsalerrors = mapsalerrors;
                        trim_cal_SAL = cal_SAL;
                        trim_cal_SAL_err = cal_SAL_err;
                        
                        jj = find(isnan(la_ptmp) == 1);
                        trimPRES(jj) = NaN;
                        trimSAL(jj) = NaN;
                        trimPTMP(jj) = NaN;
                        trim_mapped_sal(jj) = NaN;
                        trim_mapsalerrors(jj) = NaN;
                        trim_cal_SAL(jj) = NaN;
                        trim_cal_SAL_err(jj) = NaN;
                        
                        if( isempty(use_theta_lt) == 0  & isempty(use_theta_gt) == 1 )
                            jj = find(trimPTMP > use_theta_lt);
                            trimPRES(jj) = NaN;
                            trimSAL(jj) = NaN;
                            trimPTMP(jj) = NaN;
                            trim_mapped_sal(jj) = NaN;
                            trim_mapsalerrors(jj) = NaN;
                            trim_cal_SAL(jj) = NaN;
                            trim_cal_SAL_err(jj) = NaN;
                        end % if( isempty(use_theta_lt) == 0  & isempty(use_theta_gt) == 1 )
                        
                        if( isempty(use_theta_gt) == 0 & isempty(use_theta_lt) == 1 )
                            jj = find(trimPTMP < use_theta_gt);
                            trimPRES(jj) = NaN;
                            trimSAL(jj) = NaN;
                            trimPTMP(jj) = NaN;
                            trim_mapped_sal(jj) = NaN;
                            trim_mapsalerrors(jj) = NaN;
                            trim_cal_SAL(jj) = NaN;
                            trim_cal_SAL_err(jj) = NaN;
                        end % if( isempty(use_theta_gt) == 0 & isempty(use_theta_lt) == 1 )
                        
                        if( isempty(use_theta_gt) == 0 & isempty(use_theta_lt) == 0 )
                            if(use_theta_gt > use_theta_lt) 
                                % the middle band is excluded
                                jj = find(trimPTMP < use_theta_gt & trimPTMP > use_theta_lt);
                            else % if(use_theta_gt > use_theta_lt) 
                                jj = find(trimPTMP < use_theta_gt | trimPTMP > use_theta_lt);
                            end % if(use_theta_gt > use_theta_lt) 
                            trimPRES(jj) = NaN;
                            trimSAL(jj) = NaN;
                            trimPTMP(jj) = NaN;
                            trim_mapped_sal(jj) = NaN;
                            trim_mapsalerrors(jj) = NaN;
                            trim_cal_SAL(jj) = NaN;
                            trim_cal_SAL_err(jj) = NaN;
                        end % if( isempty(use_theta_gt) == 0 & isempty(use_theta_lt) == 0 )
                        
                        if( isempty(use_pres_lt) == 0 & isempty(use_pres_gt) == 1 )
                            jj = find(trimPRES > use_pres_lt);
                            trimPRES(jj) = NaN;
                            trimSAL(jj) = NaN;
                            trimPTMP(jj) = NaN;
                            trim_mapped_sal(jj) = NaN;
                            trim_mapsalerrors(jj) = NaN;
                            trim_cal_SAL(jj) = NaN;
                            trim_cal_SAL_err(jj) = NaN;
                        end % if( isempty(use_pres_lt) == 0 & isempty(use_pres_gt) == 1 )
                        
                        if( isempty(use_pres_gt) == 0 & isempty(use_pres_lt) == 1 )
                            jj = find(trimPRES < use_pres_gt);
                            trimPRES(jj) = NaN;
                            trimSAL(jj) = NaN;
                            trimPTMP(jj) = NaN;
                            trim_mapped_sal(jj) = NaN;
                            trim_mapsalerrors(jj) = NaN;
                            trim_cal_SAL(jj) = NaN;
                            trim_cal_SAL_err(jj) = NaN;
                        end % if( isempty(use_pres_gt) == 0 & isempty(use_pres_lt) == 1 )
                        
                        if( isempty(use_pres_gt) == 0 & isempty(use_pres_lt) == 0 )
                            if(use_pres_gt > use_pres_lt)
                                %  the middle band is excluded
                                jj = find(trimPRES < use_pres_gt & trimPRES > use_pres_lt);
                            else % if(use_pres_gt > use_pres_lt)
                                jj = find(trimPRES < use_pres_gt | trimPRES > use_pres_lt);
                            end % if(use_pres_gt > use_pres_lt)
                            trimPRES(jj) = NaN;
                            trimSAL(jj) = NaN;
                            trimPTMP(jj) = NaN;
                            trim_mapped_sal(jj) = NaN;
                            trim_mapsalerrors(jj) = NaN;
                            trim_cal_SAL(jj) = NaN;
                            trim_cal_SAL_err(jj) = NaN;
                        end % if( isempty(use_pres_gt) == 0 & isempty(use_pres_lt) == 0 )
                        
                        for i = 1:n
                            for j = tplot
                                if(tlevels(j) < max(trimPTMP(:, i)) ...
                                        & tlevels(j) > min(trimPTMP(:, i)))
                                    diffTheta = abs(trimPTMP(:, i) - tlevels(j));
                                    if isempty(find(~isnan(diffTheta)))
                                        Thetalevel_indexes(j, i) = NaN;
                                    else % if isempty(find(~isnan(diffTheta)))
                                        Thetalevel_indexes(j, i) = min(find(diffTheta ...
                                            == min(diffTheta)));
                                    end % if isempty(find(~isnan(diffTheta)))
                                end % if(tlevels(j) < max(trimPTMP(:, i)) ...
                            end % for j = tplot
                        end % for i = 1:n
                        
                        for i = tplot 
                            % build the S matrix for plotting
                            for j = 1:n
                                ti = Thetalevel_indexes(i, j);
                                if ~isnan(ti)
                                    % interval is one above and one below ti
                                    interval = max(ti-1,1):min(ti+1,m);
                                    a = trimPTMP(ti, j) - trimPTMP(interval, j);
                                    if( trimPTMP(ti, j) > tlevels(i) )
                                        gg = find(a > 0);
                                        if( ~isempty(gg) )
                                            % find the level with min +ve diff
                                            b = find(a == min(a(gg)));
                                            ki = interval(b);
                                        else % if( ~isempty(gg) )
                                            ki = ti;
                                        end % if( ~isempty(gg) )
                                    end % if( trimPTMP(ti, j) > tlevels(i) )
                                    if( trimPTMP(ti, j) < tlevels(i) )
                                        gg = find(a < 0);
                                        if( ~isempty(gg) )
                                            % find the level with min -ve diff
                                            b = find(-a == min(-a(gg)));
                                            ki = interval(b);
                                        else % if( ~isempty(gg) )
                                            ki=ti;
                                        end % if( ~isempty(gg) )
                                    end % if( trimPTMP(ti, j) < tlevels(i) )
                                    if( trimPTMP(ti, j) == tlevels(i) )
                                        ki = ti;
                                    end % if( trimPTMP(ti, j) == tlevels(i) )
                                    if( ki ~= ti ...
                                            & ~isnan(trimSAL(ti, j)) ...
                                            & ~isnan(trimSAL(ki, j)) ...
                                            & ~isnan(trimPTMP(ti, j)) ...
                                            & ~isnan(trimPTMP(ki, j)) )
                                        Sint(i, j) = interp1([trimPTMP(ti, j), trimPTMP(ki, j)], ...
                                            [trimSAL(ti, j), trimSAL(ki, j)], tlevels(i) );
                                    else % if( ki ~= ti ...
                                        % interpolate if possible because that is more accurate than using closest points
                                        Sint(i, j) = trimSAL(ti, j);
                                    end % if( ki ~= ti ...
                                    if( ki ~= ti & ~isnan(trim_mapped_sal(ti, j)) ...
                                            & ~isnan(trim_mapped_sal(ki ,j)) ...
                                            & ~isnan(trimPTMP(ti ,j)) ...
                                            & ~isnan(trimPTMP(ki, j)))
                                        Smap(i, j) = interp1([trimPTMP(ti, j), trimPTMP(ki, j)], ...
                                            [trim_mapped_sal(ti, j), trim_mapped_sal(ki, j)], ...
                                            tlevels(i) );
                                        Smaperr(i, j) = interp1([trimPTMP(ti, j), trimPTMP(ki, j)], ...
                                            [trim_mapsalerrors(ti, j), trim_mapsalerrors(ki, j)], ...
                                            tlevels(i) );
                                    else % if( ki ~= ti & ~isnan(trim_mapped_sal(ti, j)) ...
                                        % interpolate if possible because that is more accurate than using closest points
                                        Smap(i, j) = trim_mapped_sal(ti, j); 
                                        Smaperr(i, j) = trim_mapsalerrors(ti, j);
                                    end % if( ki ~= ti & ~isnan(trim_mapped_sal(ti, j)) ...
                                    if( ki ~= ti & ~isnan(trim_cal_SAL(ti, j)) ...
                                            & ~isnan(trim_cal_SAL(ki, j))...
                                            & ~isnan(trimPTMP(ti, j)) ...
                                            & ~isnan(trimPTMP(ki, j)) )
                                        Scal(i ,j) = interp1([trimPTMP(ti, j), trimPTMP(ki, j)], ...
                                            [trim_cal_SAL(ti, j), trim_cal_SAL(ki, j)], tlevels(i));
                                        Scalerr(i, j) = interp1([trimPTMP(ti, j), trimPTMP(ki, j)], ...
                                            [trim_cal_SAL_err(ti ,j), trim_cal_SAL_err(ki, j)], ...
                                            tlevels(i) );
                                    else % if( ki ~= ti & ~isnan(trim_cal_SAL(ti, j)) ...
                                        % interpolate if possible because that is more accurate than using closest points
                                        Scal(i, j) = trim_cal_SAL(ti, j); 
                                        Scalerr(i, j) = trim_cal_SAL_err(ti, j); 
                                    end % if( ki ~= ti & ~isnan(trim_cal_SAL(ti, j)) ...
                                end % if ~isnan(ti)
                            end % for j = 1:n
                        end % for i = tplot 
                        
% [DS2025/]
%                        figure;
                        if(goHeadless)
                            figure('Visible', 'off');
                        else % if(goHeadless)
                            figure('Visible', 'on');
                        end % if(goHeadless)
% [/DS2025]
                        set(gcf, 'defaultaxeslinewidth', 2)
                        set(gcf, 'defaultlinelinewidth', 2)
                        set(gcf, 'defaultaxesfontsize', 16)
                        
                        for k = 1:p
                            j = tplot(k);
                            subplot(p, 1, k)
                            if(isempty(Sint) == 0)
                                plot(PROFILE_NO, Sint(j, :), 'b*-');
                                hold on
                                plot(PROFILE_NO,Smap(j, :), 'r');
                                plot(PROFILE_NO,Scal(j, :), 'g');
                                mm = find(isfinite(Scal(j, :)) == 1);
                                ll = PROFILE_NO;
                                ll = ll(mm);
                                kk = Scal(j, mm);
                                nn = Scalerr(j, mm);
                                h = fill([ll, fliplr(ll)], [kk + nn, fliplr([kk - nn])], 'g');
                                set(h, 'EdgeColor', 'g');
                                errorbar(PROFILE_NO,Smap(j, :), Smaperr(j, :), 'r-')
                                plot(PROFILE_NO,Sint(j, :), 'b*-');
                                SMIN = min([Sint(j, :), Scal(j, :), Smap(j, :)]);
                                SMAX = max([Sint(j, :), Scal(j, :), Smap(j, :)]);
                                if isfinite(SMIN) & isfinite(SMAX)
                                    axis([0, max(PROFILE_NO) + 1, SMIN - 0.05, SMAX + 0.05 ])
                                end % if isfinite(SMIN) & isfinite(SMAX)
                                set(gca, 'FontSize', 12)
                                title( strcat(title_floatname, ...
                                    ' salinities with error on \theta= ', num2str(tlevels(j)), ...
                                    '^{\circ}C' ) );
                                ylabel('PSS-78')
                            end % if(isempty(Sint) == 0)
                        end % for k = 1:p
                        legend('uncal float', 'mapped salinity', 'cal float w/1xerr.', ...
                            'Location', 'Best');
                        xlabel('float profile number');
                        
                        drawnow
% [DS2025/]
% - matlab indicates that 'paperorientation' is no longer a valid parameter
%                        set(gcf, 'papertype', 'usletter', 'paperunits', 'inches', ...
%                            'paperorientation', 'portrait', 'paperposition', [0.25, 0.75, 8.0, 9.5]);
                        set(gcf, 'papertype', 'usletter', 'paperunits', 'inches', ...
                            'paperposition', [0.25, 0.75, 8.0, 9.5]);
% [/DS2025]

% [DS2025/]
%                        % CC changes 06/23 figure 6 
%                        if length(unique_cal(unique_cal > 0)) == 1
%                            print('-depsc', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, ...
%                                pn_float_dir, pn_float_name, '_6.eps'));
%                        else % if length(unique_cal(unique_cal > 0)) == 1
%                            print('-depsc', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, ...
%                                pn_float_dir, pn_float_name, '_6_split_', num2str(unique_cal(iy)), ...
%                                '.eps'));
%                        end % if length(unique_cal(unique_cal > 0)) == 1
	                    hasCtd = strfind(...
                            po_system_configuration.CONFIG_WMO_BOXES, 'ctd');
                        if(~isempty(hasCtd))
		                    if length(unique_cal(unique_cal > 0)) == 1
			                    tailtxt = '_6_ctd';
		                    else % if length(unique_cal(unique_cal > 0)) == 1
			                    tailtxt = '_6_split_ctd';
                            end % if length(unique_cal(unique_cal > 0)) == 1
                            fname = strcat(...
                                po_system_configuration.FLOAT_PLOTS_DIRECTORY, ...
                                pn_float_dir, pn_float_name, tailtxt);
                        else % if(~isempty(hasCtd))
		                    if length(unique_cal(unique_cal > 0)) == 1
                                tailtxt = '_6';
                            else % if length(unique_cal(unique_cal > 0)) == 1
                                tailtxt = '_6_split';
                            end % if length(unique_cal(unique_cal > 0)) == 1
                            fname = strcat(...
                                po_system_configuration.FLOAT_PLOTS_DIRECTORY, ...
                                pn_float_dir, pn_float_name, tailtxt);
                        end % if(~isempty(hasCtd))
                        print(pltFileType, '-r300', fname)
% [/DS2025]
                    end % if unique_cal(iy) > 0
                end % for iy = 1:n_seq

                % Brian King's plot: salinity anomaly time series on theta levels (figure 7) ------------
% [DS2025/]
%                figure;
                if(goHeadless)
                    figure('Visible', 'off');
                else
                    figure('Visible', 'on');
                end % if(goHeadless)
% [/DS2025]
                set(gcf, 'defaultaxeslinewidth', 2)
                set(gcf, 'defaultlinelinewidth', 2)
                set(gcf, 'defaultaxesfontsize', 16)

                fl.useqc = '0';
                fl.plot = 1;
                fl.yaxes = [2 5 20];
                d.PSAL = cal_SAL;
                d.TEMP = TEMP;
                d.PRES = PRES;
                d.PSAL_QC = zeros(m, n);
                d.TEMP_QC = zeros(m, n);
                d.PRES_QC = zeros(m, n);
                d.LONGITUDE = LONG;
                d.LATITUDE = LAT;
                d.PROFILE_NO = PROFILE_NO;
                fl = anom(d, fl); % Brian King's routine
                subplot('position', [0.1 0.45 0.8 0.35])
                title(['Calibrated salinity anom on theta. ' title_floatname])

                drawnow
% [DS2025/]
% - matlab indicates that 'paperorientation' is no longer a valid parameter
%                set(gcf, 'papertype', 'usletter', 'paperunits', 'inches',' paperorientation',...
%                    'portrait', 'paperposition', [0.25, 0.5, 8.0, 10.0]);
                set(gcf, 'papertype', 'usletter', 'paperunits', 'inches', ...
                    'paperposition', [0.25, 0.5, 8.0, 10.0]);
% [/DS2025]

% [DS2025/]
%                print('-depsc', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, ...
%                    pn_float_dir, pn_float_name, '_7.eps'));
                hasCtd = strfind(po_system_configuration.CONFIG_WMO_BOXES, ...
                    'ctd');
                if(~isempty(hasCtd))
                    fname = strcat(...
                        po_system_configuration.FLOAT_PLOTS_DIRECTORY, ...
                        pn_float_dir, pn_float_name, '_7_ctd');
                else % if(~isempty(hasCtd))
                    fname = strcat(...
                        po_system_configuration.FLOAT_PLOTS_DIRECTORY, ...
                        pn_float_dir, pn_float_name, '_7');
                end % if(~isempty(hasCtd))
                print(pltFileType, '-r300', fname);
% [/DS2025]

                % Paul Robbins' analyse variance plot (figure 8) ------------
                % CC changes 06/23 figure 8  to take into account the splitting of the time series 
                % add l647 to l666
                unique_cal = unique(lo_float_calseries.calseries);
                n_seq = length(unique_cal);

                for iy = 1:n_seq
                    calindex = find(lo_float_calseries.calseries == unique_cal(iy));
                    k = length(calindex);
                    if unique_cal(iy) > 0
                        % choose 10 float theta levels to use in the piecewise linear fit --------
                        unique_SAL = SAL(:, calindex);
                        unique_cal_SAL = cal_SAL(:, calindex);
                        unique_PTMP = PTMP(:, calindex);
                        unique_PRES = PRES(:, calindex);
                        unique_la_ptmp = la_ptmp(:, calindex);
                        unique_mapped_sal = mapped_sal(:, calindex);
                        unique_mapsalerrors = mapsalerrors(:, calindex);
                        
                        
                        [tlevels, plevels, index, var_s_Thetalevels, Thetalevels] = find_10thetas(...
                            unique_SAL, unique_PTMP, unique_PRES, unique_la_ptmp, use_theta_gt, ...
                            use_theta_lt, use_pres_gt, use_pres_lt, use_percent_gt);  
                        
% [DS2025/]
%                        figure;
                        if(goHeadless)
                            figure('Visible', 'off');
                        else
                            figure('Visible', 'on');
                        end % if(goHeadless)
% [/DS2025]
                        set(gcf, 'defaultaxeslinewidth', 1)
                        set(gcf, 'defaultlinelinewidth', 1)
                        set(gcf, 'defaultaxesfontsize', 12)
                        
                        % plot t-s profile
                        subplot(222)
                        %plot(SAL,PTMP,'b-');
                        plot(unique_SAL,unique_PTMP, 'b-');  %CC changes 06/23
                        ylabel('Potential temp (^{\circ}C)')
                        x = get(gca, 'xlim');
                        y = get(gca, 'ylim');
                        xlabel('PSS-78')
% [DS2025/]
% - it apprears strcat() strips trailing blanks
%                        title(strcat('OW chosen levels - ', pn_float_name));
                        title(['OW chosen levels - ', pn_float_name]);
% [/DS2025]
                        for i = 1:10
                            hold on
                            plot(x, [tlevels(i) tlevels(i)], 'g-');
                        end % for i = 1:10
                        
                        % plot s var on t
                        subplot(221)
                        plot(var_s_Thetalevels, Thetalevels, 'b-')
                        x = get(gca, 'xlim');
                        for i = 1:10
                            hold on
                            plot(x, [tlevels(i) tlevels(i)] ,'g-');
                        end % for i = 1:10
                        %title('Salinity Variance on Theta')
% [DS2025/]
%                        title(['Salinity Variance on Theta, cycles ' ...
%                            num2str(PROFILE_NO(calindex(1))) '-' ...
%                            num2str(PROFILE_NO(calindex(end))) ])  %CC changes 06/23
                        title({'Salinity Variance on Theta', ...
                            ['cycles: ' num2str(PROFILE_NO(calindex(1))) '-' ...
                            num2str(PROFILE_NO(calindex(end)))]})
% [/DS2025]
                        ylabel('Potential temp (^{\circ}C)')
                        xlabel('salinity variance')
                        %set(gca,'ylim',y)
                        
                        % plot p-t profile
                        subplot(223)
                        %plot(PTMP,-PRES,'b-'); 
                        plot(unique_PTMP, -unique_PRES, 'b-'); %CC changes 06/23
                        x = get(gca,'xlim');
                        xlabel('^{\circ}C')
                        ylabel('Pressure (dbar)')
% [DS2025/]
%                        title(strcat('OW chosen levels - ', pn_float_name));
                        title(['OW chosen levels - ', pn_float_name]);
% [/DS2025]
                        for i = 1:10
                            hold on
                            plot(x,[-plevels(i) -plevels(i)] ,'g-');
                        end % for i = 1:10
                        
                        % plot p-s profile
                        subplot(224)
                        %plot(SAL,-PRES,'b-');
                        plot(unique_SAL,-unique_PRES, 'b-'); %CC changes 06/23
                        x = get(gca, 'xlim');
                        xlabel('PSS-78')
% [DS2025/]
%                        title(strcat('OW chosen levels - ', pn_float_name));
                        title(['OW chosen levels - ', pn_float_name]);
% [/DS2025]
                        for i = 1:10
                            hold on
                            plot(x, [-plevels(i) -plevels(i)], 'g-');
                        end % for i = 1:10
                        
                        drawnow
% [DS2025/]
% - matlab indicates that 'paperorientation' is no longer a valid parameter
%                        set(gcf, 'papertype', 'usletter', 'paperunits', 'inches', ...
%                            'paperorientation', 'portrait', 'paperposition', ...
%                            [0.5, 0.25, 8.0, 10.25]);
                        set(gcf, 'papertype', 'usletter', 'paperunits', 'inches', ...
                            'paperposition', [0.5, 0.25, 8.0, 10.25]);
% [/DS2025]

% [DS2025/]
%                         if length(unique_cal(unique_cal > 0)) == 1
%                             print('-depsc', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, ...
%                                 pn_float_dir, pn_float_name, '_8.eps'));
%                         else % if length(unique_cal(unique_cal > 0)) == 1
%                             print('-depsc', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, ...
%                                 pn_float_dir, pn_float_name, '_8_split_', num2str(unique_cal(iy)), ...
%                                 '.eps'));
                        end % if length(unique_cal(unique_cal > 0)) == 1
                        hasCtd = strfind(...
                            po_system_configuration.CONFIG_WMO_BOXES, 'ctd');
                        if(~isempty(hasCtd))
		                    if length(unique_cal(unique_cal > 0)) == 1
                			    tailtxt = '_8_ctd';
                            else % if length(unique_cal(unique_cal > 0)) == 1
			                    tailtxt = '_8_split_ctd';
		                    end % if length(unique_cal(unique_cal > 0)) == 1
                            fname = strcat(...
                                po_system_configuration.FLOAT_PLOTS_DIRECTORY, ...
                                pn_float_dir, pn_float_name, tailtxt);
                        else % if(~isempty(hasCtd))
                            if length(unique_cal(unique_cal > 0)) == 1
                                tailtxt = '_8';
                            else % if length(unique_cal(unique_cal > 0)) == 1
                                tailtxt = '_8_split';
                            end % if length(unique_cal(unique_cal > 0)) == 1
                            fname = strcat(...
                                po_system_configuration.FLOAT_PLOTS_DIRECTORY, ...
                                pn_float_dir, pn_float_name, tailtxt);
                        end % if(~isempty(hasCtd))
                        print(pltFileType, '-r300', fname);
% [/DS2025]
                    end % if unique_cal(iy) > 0
                end % for iy = 1:n_seq
            end %if(isempty(find(isnan(cal_SAL) == 0)) == 0) ---------------
        end %if(isempty(find(isnan(mapped_sal) == 0)) == 0) -----------
    end %if(isempty(find(isnan(PRES) == 0)) == 0) ------------------


