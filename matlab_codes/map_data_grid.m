
function [vgrid,vgriderror,vdata,vdataerror]=map_data_grid(v,posgrid,posdata,gs,ts,q,s,e,phi,pv)

% **********************************************************************
% Objective mapping routine - needs covarxyt_pv.m.
%
% Breck Owens, November 2007
%
% Inputs:
%
% v=data   rank*m1 matrix of data (rank is the number of repeated data). This does not take NaN's.
% posgrid  m2*4 matrix of grid locations [latitude,longitude,dates,depth]
% posdata  m1*4 matrix of data locations [latitude,longitude,dates,depth]
% gs       scalar longitude scale (in degrees)
% ts       scalar latitude scale (in degrees)
% q        scalar time scale (in years)
% s        scalar signal level (in variance terms)
% e        scalar noise level (in variance terms)
% phi      scalar cross isobaric scale
% pv       include potential vorticity in the covariance or not (1=include, 0=exclude)
%
% Outputs:
%
% vgrid    m2*1 matrix of mapped fields (n is the number of repeated data)
% vgriderror  m2*1 matrix of error estimates of the mapped data
%             (assumed statistics)
% vdata    m1*1 matrix of mapped fields on original data locations
% vdataerror  m1*1 matric of error estimate of the mapped data
%             (assumed statitics)
%
% Purpose: An optimal mapping routine, taking data measured in arbitrary
%          geographic locations and mapping these data onto a more regular
%          grid. As should happen in every mapping problem the data are both
%          mapped onto the prescribed grid, and the data locations (so that
%          the mapped field can be checked with the original data to ensure
%          that the statistics are valid and consistent).
%          Before objectively mapping the data, a mean, using the
%          correlation scales to define the weights for the sum (see
%          Bretherton, etal, 1975) is removed.  The error estimate includes
%          the contributions from both the mapping and the mean
%          estimatation.
%
% Delphine Dobler (DD), August 2024: 
%            3.3 - Performance: Inside map_data_grid, comment calculation
%            of vdataerror as this is never used afterwards.
% **********************************************************************

[m,n]=size(v);
[nnn,i]=size(posgrid);
[mmm,i]=size(posdata);

% to make array multiplications work out -----

v = v';

vdata=zeros(mmm,1);
vdataerror=vdata;
vgrid=zeros(nnn,1);
vgriderror=vgrid;

% Create the data-data covariance matrix ----------
Cdd = inv(s*covarxyt_pv( posdata,posdata,gs,ts,q, phi, pv)+ e*eye(mmm,mmm));

sum_Cdd = sum(sum(Cdd)); % used to estimate the mean field

% estimate the mean field using the OI weights (eqn 21 of Bretherton, etal)
vm = sum(Cdd*v)/sum_Cdd;

wght = Cdd*( v - vm );

% calculate the objectively mapped fields on data and grid

% map to posdata -----

% data-mapping point covariance matrix --------
Cmd = s*covarxyt_pv(posdata, posdata, gs, ts, q, phi, pv);
vdata = Cmd*wght + vm;
% include error in mean (eqn 24 of Bretherton, etal)
%vdataerror = repmat(sqrt(s-diag(Cmd*Cdd*Cmd') ...
%                         +(1- sum(Cmd*Cdd,2)).^2/sum_Cdd),1,1);
% DD (2024/08-3.3): never used afterwards, save some time (Uncomment to
% compute)
vdataerror = [];  

% map to posgrid -----

Cmd = s*covarxyt_pv(posgrid, posdata, gs, ts, q, phi, pv);
vgrid = Cmd*wght + vm;
vgriderror = repmat(sqrt(s-diag(Cmd*Cdd*Cmd') ...
			   +(1- sum(Cmd*Cdd,2)).^2/sum_Cdd),1,1);

return

