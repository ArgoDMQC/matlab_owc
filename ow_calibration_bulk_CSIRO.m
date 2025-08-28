% ow_calibration_bulk_CSIRO.m
%
% [07/07/21 Dirk Slawinski
% - lets get this script localized and portable

% float areas
flt_dir = '';
% the floats to process
%wmoids = [5905405];
wmoids = [5905025, 5905414];
% do we also for CTD?
doCTD = true;
% go headless 
goHeadless = true;
%goHeadless = false;

% some configs
% the current Matlab based OWC directory
OWCdir = '/home/argo/ArgoDM/dmqc_ow/OWC_3_dev';
currentDir = pwd();

% go to OWC 
cd(OWCdir);

% set the path
owroot = OWCdir; % '/home/argo/ArgoDM/dmqc_ow/OWC_3_dev';
addpath(genpath(owroot))

% now loop over the lo_system_configuration = load_configuration( 'ow_config.txt' );floats
for wmoid = wmoids
    % get a local version of the config fata
    lo_system_configuration = load_configuration( 'ow_config.txt' );

    %flt_dir = '';
    flt_name = sprintf('%d', wmoid);
    
    fprintf('%s Working on %d\n', datestr(now), wmoid);
    %fprintf('%s\n', lo_system_configuration.CONFIG_WMO_BOXES);
    
    % the default OWC mapped name
    genericMappedName = sprintf('%s/data/float_mapped/map_%d.mat', OWCdir, wmoid);
    % CTD name we want it to have
    ctdMappedName = sprintf('%s/data/float_mapped/map_%d_ctd.mat', OWCdir, wmoid);
    % ARGO name we want it to have
    argoMappedName = sprintf('%s/data/float_mapped/map_%d_argo.mat', OWCdir, wmoid);
    
    % sometimes we want to use the historical CTD data
    if(doCTD)
        % swap argo for ctd comparison
        lo_system_configuration.CONFIG_WMO_BOXES = ...
            strrep(lo_system_configuration.CONFIG_WMO_BOXES, 'argo', 'ctd');
        % check on mapped file status
        if(exist(ctdMappedName, 'file'))
            % we have run CTD before so lets copy it to the working/generic
            % name - copyfile(src, dest)
            % but 1st check if _argo exists.
            if(exist(argoMappedName, 'file'))
                % yes, so no need to backup
                fprintf('keeping existing _argo file during ctd comparison\n');
            elseif (exist(genericMappedName, 'file')) % if(exist(argoMappedName, 'file'))
                % rename generic as _argo ie backup
                mvStat = movefile(genericMappedName, argoMappedName);
                % if a move error happens go back
                if(~mvStat)
                    fprintf('error whilst renaming generic to %s\n', argoMappedName);
                    % now return to the previous dir
                    cd(currentDir);
                    return
                end % if(~mvStat)
                fprintf('backing up generic file as %s during ctd comparison\n', argoMappedName);
            end % if(exist(argoMappedName, 'file'))
            % ow copy the _ctd file to generic
            mvStat = copyfile(ctdMappedName, genericMappedName);
            % if a move error happens go back
            if(~mvStat)
                fprintf('error whilst copying _ctd to %s\n', genericMappedName);
                % now return to the previous dir
                cd(currentDir);
                return
            end % if(~mvStat)
            fprintf('using %s as generic file\n', ctdMappedName);
        else % if(exist(ctdMappedName, 'file'))
            % ctd does not yet exist, so assume generic file is the _argo one
            % which needs to be backed up to _argo
            % 1st check if it exists
            if(exist(genericMappedName, 'file'))
                mvStat = movefile(genericMappedName, argoMappedName);
                % if a move error happens go back
                if(~mvStat)
                    fprintf('error whilst renaming generic to %s\n', argoMappedName);
                    % now return to the previous dir
                    cd(currentDir);
                    return
                end % if(~mvStat)
                fprintf('backing up generic file as %s during ctd comparison\n', argoMappedName);
            else % if(exist(genericMappedName, 'file'))
                % does not exist let it create the generic
                fprintf('no existing generic file during ctd comparison\n');
            end % if(exist(genericMappedName, 'file'))
        end % if(exist(ctdMappedName, 'file'))
        % at this point the generic mapped name should either be the _ctd one or
        % none at all, 1st ctd run
        
        % now run it again
        update_salinity_mapping( flt_dir, flt_name, lo_system_configuration );
        set_calseries( flt_dir, flt_name, lo_system_configuration );
        calculate_piecewisefit( flt_dir, flt_name, lo_system_configuration );
        plot_diagnostics_ow_CSIRO( flt_dir, flt_name, lo_system_configuration, goHeadless );
        
        % swap ctd back to argo
        lo_system_configuration.CONFIG_WMO_BOXES = ...
            strrep(lo_system_configuration.CONFIG_WMO_BOXES, 'ctd', 'argo');
        
        % backup the _ctd mapped data 
        mvStat = movefile(genericMappedName, ctdMappedName);
        % if a move error happens go back
        if(~mvStat)
            fprintf('error whilst renaming to %s\n', ctdMappedName);
            % now return to the previous dir
            cd(currentDir);
            return
        end % if(~mvStat)
        fprintf('generic file saved as %s during ctd comparison\n', ctdMappedName);
    end % if(doCTD)
    
    % let's use the ARGO boxes data
    % swap ctd to argo, just to be safe
    lo_system_configuration.CONFIG_WMO_BOXES = ...
        strrep(lo_system_configuration.CONFIG_WMO_BOXES, 'ctd', 'argo');
    
    % check on mapped file status
    if(exist(argoMappedName, 'file'))
        % we have an _argo, copy it to the working/generic name - 
        % copyfile(src, dest)
        mvStat = copyfile(argoMappedName, genericMappedName);
        % if a move error happens go back
        if(~mvStat)
            fprintf('error whilst copying _argo to %s\n', genericMappedName);
            % now return to the previous dir
            cd(currentDir);
            return
        end % if(~mvStat)
        fprintf('using %s as generic file during argo comparison\n', argoMappedName);
    else % if(exist(argoMappedName, 'file'))
        % does not exist so the generic is likely the _argo data
        % NOTE: it could be the generic or could be the _ctd but let's be
        % brave
        fprintf('using existing generic file during argo comparison\n');
    end % if(exist(argoMappedName, 'file'))

    update_salinity_mapping( flt_dir, flt_name, lo_system_configuration );
    set_calseries( flt_dir, flt_name, lo_system_configuration );
    calculate_piecewisefit( flt_dir, flt_name, lo_system_configuration );
    plot_diagnostics_ow_CSIRO( flt_dir, flt_name, lo_system_configuration, goHeadless );

    % backup the _argo mapped data and keep generic as _argo
    mvStat = copyfile(genericMappedName, argoMappedName);
    % if a move error happens go back
    if(~mvStat)
        fprintf('error whilst copying to %s\n', argoMappedName);
        % now return to the previous dir
        cd(currentDir);
        return
    end % if(~mvStat)
    fprintf('generic file saved as %s during argo comparison\n', argoMappedName);
end % wmoid

% now return to the previous dir
fprintf('All Done!\n');
cd(currentDir);
return
