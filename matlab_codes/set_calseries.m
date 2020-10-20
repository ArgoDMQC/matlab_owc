
function set_calseries( pn_float_dir, pn_float_name, po_system_configuration )

%
% function set_calseries( pn_float_dir, pn_float_name, po_system_configuration )
%
% Annie Wong, September 2008
% Breck Owens, October 2006
%


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
    display(['***tip: to modify these parameters, first delete the file: ' ls_calseries_filename ' and edit set_calseries.m'])
catch
    breaks = [];
    max_breaks = 3;   % 0 for linear trend and -1 for offset only!
    calseries = [ones(1,n)];
    % calseries = [ones(1,33) 2*ones(1,n-33)];    % example: split the time series at profile 33
    % calseries = [ones(1,33) 0  ones(1,n-33-1)]; % example: ignore profile 34
    calib_profile_no = PROFILE_NO;
    use_theta_lt = [];
    use_theta_gt = [];
    use_pres_gt = [];
    use_pres_lt = [];
    use_percent_gt = 0.5;
end

display(['calseries: ' num2str(calseries)])
display(['breaks = ' num2str(breaks)])
display(['max_breaks = ' num2str(max_breaks)])
display(['use_theta_lt = ' num2str(use_theta_lt)])
display(['use_theta_gt = ' num2str(use_theta_gt)])
display(['use_pres_gt = ' num2str(use_pres_gt)])
display(['use_pres_lt = ' num2str(use_pres_lt)])
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

cn = length(calib_profile_no);

for i=1:length(missing_profile_index)
   j = missing_profile_index(i);
   calib_profile_no = [calib_profile_no, PROFILE_NO(j)];
   calseries = [calseries, calseries(max(j-1,1))]; % same flag as previous profile
end


% sort the calseries file by profile_number ----

[y,ii]=sort(calib_profile_no);

calib_profile_no=calib_profile_no(ii);
calseries=calseries(ii);


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
  end
end


% save calseries file ----

save( ls_calseries_filename, 'breaks', 'max_breaks', 'calseries', 'calib_profile_no', 'use_theta_lt', 'use_theta_gt', 'use_pres_gt', 'use_pres_lt', 'use_percent_gt' );

