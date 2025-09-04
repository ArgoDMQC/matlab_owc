# plot\_diagnostics\_ow\_CSIRO(...)

## Introduction

This version of *plot\_diagnostics\_ow* has been modified by Dirk Slawinski from *CSIRO* to work with newer Matlab versions. Starting with the 2018 version of Matlab many handle graphics behaviors changed the way plotting happened particularly for grouped plot items. Each plot item is now its own entity and so this requires rewrites to ensure plots look and behave the same as in previous versions of Matlab.

Using the new routine with the previous parameters works as before, plots no longer have large legends nor require 3D rendering to achieve x-error bars in Figs 2 and 4 and maintain their interactivity. 

Brian King's *anom.m* has also been undated so you no longer see the ghosts of the plotted over xlabels from the top half of the main plot.
This version also incorporates Delphine Dobler's splitting code.
Plots are saved as PNGs by default [1].

	plot_diagnostics_ow_CSIRO( pn_float_dir, pn_float_name, po_system_configuration)
	
A new feature added to the *CSIRO* variant is passing in the parameter 'goHeadless' at the end.  It is a logical, TRUE/FALSE, with default FALSE if not provided. This forces rendering straight to files rather than providing interactive figure plots.

	plot_diagnostics_ow_CSIRO( pn_float_dir, pn_float_name, po_system_configuration, goHeadless)


[1] this can be changed in the script by setting:

	pltFileType = '-dpng';

to

	pltFileType = '-depsc';

One could potentially also use '-djpeg' or even '-dsvg' but PNGs are better overall for figures and SVGs may be cumbersome for post editing.

## Usage
### Traditional
The traditional way of calling this routine is as the last call in the *ow\_calibration.m* script:

	for i=1:length(float_names)

    	flt_dir = float_dirs{i};
    	flt_dir = deblank(flt_dir);
    	flt_name = float_names{i};

    	disp([datestr(now) ' Working on ' flt_name])

    	update_salinity_mapping( flt_dir, flt_name, lo_system_configuration );

    	set_calseries( flt_dir, flt_name, lo_system_configuration );

    	calculate_piecewisefit( flt_dir, flt_name, lo_system_configuration );

    	plot_diagnostics_ow( flt_dir, flt_name, lo_system_configuration );
	end % for i=1:length(float_names)

	
To use the *CSIRO* verion just replace *plot\_diagnostics\_ow* with *plot\_diagnostics\_ow\_CSIRO* and all should be fine.

	plot_diagnostics_ow_CSIRO( flt_dir, flt_name, lo_system_configuration );


### Headless
Passing in the new added parameter *goHeadless* as logical true (or 1) the routine will not render to screen but straight to image/vector files. For file types see [1] in **Introduction**. This means no interactive plots which is useful in bulk processing. If the parameter is missing, the routine renders interactive plots. If it is present but false or 0 then the routine renders interactive plots.

### Bulk Processing
The original calibration script supports bulk processing, but it's limited by how many floats can be reviewed at once. To address this, *CSIRO* introduced a headless mode, which also simplifies automated documentation generation.

*CSIRO* has also provided an enhanced script, *ow\_calibration\_bulk\_CSIRO.m*, which processes multiple floats like the original but adds support for calibrating against both mapped CTD data and historical ARGO data. Configuration remains at the top of the script, so the usage is largely unchanged.


	% float areas
	flt_dir = '';
	% the floats to process
	wmoids = [5905025, 5905414];
	% do we also for CTD?
	doCTD = true;
	% go headless 
	goHeadless = true;
	%goHeadless = false;








