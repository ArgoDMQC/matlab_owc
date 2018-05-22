
%
% have to edit "ow_config.txt" if this system is moved somewhere else.
% ----------------------------------------------------------------------
%
% script ow_calibration
%
% dir names and float names have to correspond,
% e.g. float_dirs={'sio/';'uw/'}
%      float_names={'R49000139';'R39033'};
%
% these variables have to be set before ow_calibration is called.
%

lo_system_configuration = load_configuration( 'ow_config.txt' );

for i=1:length(float_names)

    flt_dir = float_dirs{i};
    flt_dir = deblank(flt_dir);
    flt_name = float_names{i};

  disp([datestr(now) ' Working on ' flt_name])

  update_salinity_mapping( flt_dir, flt_name, lo_system_configuration );

  set_calseries( flt_dir, flt_name, lo_system_configuration );

  calculate_piecewisefit( flt_dir, flt_name, lo_system_configuration );

  plot_diagnostics_ow( flt_dir, flt_name, lo_system_configuration );

end

