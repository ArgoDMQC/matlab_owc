
function set_calseries( pn_float_dir, pn_float_name, po_system_configuration )

%
% function set_calseries( pn_float_dir, pn_float_name, po_system_configuration )
%
% Annie Wong, September 2008
% Breck Owens, October 2006
% Delphine Dobler (DD), August 2024: 2 - create output directory if they do not exist already
% Delphine Dobler (DD), June 2026: 
%       - separate continuous piecewisefit and multiple theta sets segmentation
%       - expand comments for the configuration parameters
%       - remove unused parameter cn


% DD (2024/08-2)
outdir=[ po_system_configuration.FLOAT_CALIB_DIRECTORY pn_float_dir ];
if not(isfolder(outdir))
    mkdir(outdir)
end

% load data ---

lo_float_source_data = load( strcat( po_system_configuration.FLOAT_SOURCE_DIRECTORY, pn_float_dir, pn_float_name, po_system_configuration.FLOAT_SOURCE_POSTFIX ) ) ;

PROFILE_NO = lo_float_source_data.PROFILE_NO;
n=length(PROFILE_NO);

ls_calseries_filename = strcat( po_system_configuration.FLOAT_CALIB_DIRECTORY, pn_float_dir, po_system_configuration.FLOAT_CALSERIES_PREFIX, pn_float_name, po_system_configuration.FLOAT_CALIB_POSTFIX ) ;


% build default values ---
disp(' ')
disp('___________________________________________')
disp('SET CALSERIES PARAMETERS')

try
    load(ls_calseries_filename);
    disp(['***tip: to modify these parameters, first delete the file: ' ls_calseries_filename ' and edit set_calseries.m'])
catch
    
    % 1a - prescribe the calseries for the theta levels selection (e.g. front crossing)
    calseries = [ones(1,n)];
    %calseries = [ones(1,56) 2*ones(1,n-56)]; % example: two theta sets splitting at cycle 56
    % calseries = [ones(1,33) 0  ones(1,n-33-1)]; % example: ignore profile 34
    
    % 1b - further configure the theta levels selection
    use_theta_lt = [];
    use_theta_gt = [];
    use_pres_gt = [];
    use_pres_lt = [];
    use_percent_gt = 0.5;
    
    % 2 - prescribe parameterisation for the piecewisefit
    
    % 2a - First set the continuous segments (no jump)
    calseries_for_cont_pwf = calseries;
    %calseries_for_cont_pwf  = [ones(1,n)]; % example: even though two theta sets are defined, ensure the proposed adjustment is continuous
    
    % 2b - prescribe some manual breaks
    breaks = []; % don't bother anymore when several segments: list all the breaks you wish (e.g. [30, 40, 120])
    
    % 2c - prescribe the maximum number of breaks by piecewisefit segment
    % Either a scalar or an array (e.g. [1,2,2] if 3 segments in calseries_for_cont_pwf). 
    % If a scalar or an array of dimension 1, the same maximal number of breaks is applied to each segment.
    % If an array of dimension > 1, the array dimension should match the number of segments in calseries_for_cont_pwf
    max_breaks = [];   % 0 for linear trend and -1 for offset only.
    
    
    calib_profile_no = PROFILE_NO;
    
end

disp(['calseries (theta levels selection): ' num2str(calseries)])
disp(['calseries (piecewisefit with continuous adj sections): ' num2str(calseries_for_cont_pwf)])
disp(['breaks = ' num2str(breaks)])
disp(['max_breaks = ' num2str(max_breaks)])
disp(['use_theta_lt = ' num2str(use_theta_lt)])
disp(['use_theta_gt = ' num2str(use_theta_gt)])
disp(['use_pres_gt = ' num2str(use_pres_gt)])
disp(['use_pres_lt = ' num2str(use_pres_lt)])
disp('___________________________________________')
disp(' ')

% to enhance backward compatiability because I added a new variable "use_percent_gt" and changed 99999 to [] in Sep08 ---

if(exist('use_percent_gt')==0)
  use_percent_gt = 0.5;
end

if use_theta_gt == 99999; use_theta_gt = [];  end
if use_theta_lt == 99999; use_theta_lt = [];  end
if use_pres_gt == 99999; use_pres_gt = [];  end
if use_pres_lt == 99999; use_pres_lt = [];  end


% compare profile_number in source file and calseries file ----

missing_profile_index = [];

for i=1:n
   a=find( calib_profile_no==PROFILE_NO(i) );
   if( isempty(a)==1 )
     missing_profile_index = [ missing_profile_index, i ];
   end
end


% update calseries by missing_profile_index ----

for i=1:length(missing_profile_index)
   j = missing_profile_index(i);
   calib_profile_no = [calib_profile_no, PROFILE_NO(j)];
   calseries = [calseries, calseries(max(j-1,1))]; % same flag as previous profile
   calseries_for_cont_pwf = [calseries_for_cont_pwf, calseries_for_cont_pwf(max(j-1,1))];
end


% sort the calseries file by profile_number ----

[~,ii]=sort(calib_profile_no);

calib_profile_no=calib_profile_no(ii);
calseries=calseries(ii);
calseries_for_cont_pwf=calseries_for_cont_pwf(ii);


% if SAL or TEMP or PRES = all NaNs, calseries = 0 -----

SAL = lo_float_source_data.SAL;
TEMP = lo_float_source_data.TEMP;
PRES = lo_float_source_data.PRES;

for i=1:n
  ii=find(isnan(SAL(:,i))==0);
  jj=find(isnan(TEMP(:,i))==0);
  kk=find(isnan(PRES(:,i))==0);
  if(isempty(ii)==1|isempty(jj)==1|isempty(kk)==1)
     calseries(i)=0;
     calseries_for_cont_pwf(i)=0;
  end
end


% save calseries file ----

save(ls_calseries_filename, 'breaks', 'max_breaks', 'calseries','calseries_for_cont_pwf', 'calib_profile_no', 'use_theta_lt', 'use_theta_gt', 'use_pres_gt', 'use_pres_lt', 'use_percent_gt' );

