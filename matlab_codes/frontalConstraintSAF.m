
function [ hist_sal, hist_ptmp, hist_pres, hist_lat, hist_lon, hist_dates, hist_Z ] = frontalConstraintSAF( grid_sal, grid_ptmp, grid_pres, grid_lat, grid_long, grid_dates, grid_Z, LAT, LONG,Z,PTMP,SAL, po_config_data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%													%
%	frontalConstraintSAF.m 										%
%													%
%	Input :												%
%		grid_sal, grid_ptmp, grid_pres, grid_lat, grid_long, grid_dates, Grid_Z: 		%
%		-> Array of the first selection of historical data					%
%		LAT,LONG,Z,PTMP,SAL: 									%
%		-> Argo profile data									%
%													%
%	Output:												%
%		hist_sal, hist_ptmp, hist_pres, hist_lat, hist_lon, hist_dates, hist_Z:			%
%		-> Array of historical data after application of the SAF frontal criterion		%
%													%
%													%
%	Created by J.B. Sallée; 15/06/2006								%
%	sallee@legos.cnes.fr										%
%													%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%grid_sal=la_bhist_sal;
%grid_ptmp=la_bhist_ptmp;
%grid_pres=la_bhist_pres;
%grid_lat=la_bhist_lat;
%grid_long=la_bhist_long;
%grid_dates=la_bhist_dates;
%grid_Z=la_bhist_Z;
%Z=PRES;
%po_config_data=po_system_configuration;

load(strcat( po_config_data.CONFIG_DIRECTORY, po_config_data.CONFIG_SAF ));
%keep only below 100m depth data :
S_meanS=S_meanS(3:end);
S_meanN=S_meanN(3:end);
S_stdS=S_stdS(3:end);
T_stdS=T_stdS(3:end);
S_stdN=S_stdN(3:end);
T_stdN=T_stdN(3:end);
T_meanS=T_meanS(3:end);
T_meanN=T_meanN(3:end);
Deph=Deph(3:end);



if (LAT<-30)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 			ARGO Float			 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
TEMP=sw_temp(SAL,PTMP,Z,0);
isok=find(~isnan(Z) &~isnan(TEMP));
if (isempty(find(diff(Z(isok))==0)) & length(isok)>2 & min(Z(isok))<300 & max(Z(isok))>300 )
	
   % SAF :
   T300=interp1(Z(isok),TEMP(isok),300);
   if ~isempty(T300)
   	if (T300>5)
		CritSAF=1;
   	elseif T300<3
		CritSAF=-1;
   	else
		CritSAF=0;
   	end
   else
	CritSAF=0;
   end

% Second step : the envelope test ...
if (CritSAF==0 )
	SAL_int=interp1(Z(isok),SAL(isok),Deph');
	TEMP_int=interp1(Z(isok),TEMP(isok),Deph');
	Southinf_int=NaN*ones(size(Deph'));
	Southsup_int=NaN*ones(size(Deph'));
	Northinf_int=NaN*ones(size(Deph'));
	Northsup_int=NaN*ones(size(Deph'));

	Northinf_int=interp1(T_meanN-T_stdN,S_meanN-S_stdN,TEMP_int);
	Northsup_int=interp1(T_meanN+T_stdN,S_meanN+S_stdN,TEMP_int);
	Southinf_int=interp1(S_meanS-S_stdS,T_meanS-T_stdS,SAL_int);
	Southsup_int=interp1(S_meanS+S_stdS,T_meanS+T_stdS,SAL_int);

	isok2=find(~isnan(SAL_int)&~isnan(Southinf_int) &~isnan(Southsup_int) &~isnan(Northinf_int) &~isnan(Northsup_int) & Deph'>150 & Deph'<1700);

	if ~isempty(isok2)	
		PtSouth=find(TEMP_int(isok2)>Southinf_int(isok2) & ...
				TEMP_int(isok2)<Southsup_int(isok2));
		PtNorth=find(SAL_int(isok2)>Northinf_int(isok2) & ...
				SAL_int(isok2)<Northsup_int(isok2));		


	isSouth=0; isNorth=0;
	if (length(PtSouth)==length(isok2))
		isSouth=1;
	end
	if (length(PtNorth)==length(isok2))
		isNorth=1;
	end	

	if (isSouth & isNorth)
		CritSAF=0;
	elseif (isSouth)
		CritSAF=-1;
	elseif (isNorth)
		CritSAF=1;
	else
		CritSAF=0;
	end
	end % if ~isempty(isok2)	
end % if (CritSAF==0 )

else % if (isempty(find(diff(Z(isok))==0)) & length(isok)>2)
	CritSAF=0;
end % if (isempty(find(diff(Z(isok))==0)) & length(isok)>2)


%%%%%%%%% Debug lines ....
Check=0;
if (Check==1)
	figure(1)
	clf
	hold on
	disp(['CritSAF : ', int2str(CritSAF)])
	plot(S_meanS,T_meanS,'r')
	plot(S_meanN,T_meanN,'k')
	plot(S_meanS+S_stdS,T_meanS+T_stdS,'r--')
	plot(S_meanS-S_stdS,T_meanS-T_stdS,'r--')
	plot(S_meanN+S_stdN,T_meanN+T_stdN,'k--')
	plot(S_meanN-S_stdN,T_meanN-T_stdN,'k--')
	plot(SAL,TEMP,'b')
	%plot(SAL_int,TEMP_int,'g--')	
	figure(2)
	plot(SAL_int,TEMP_int,'b')
	plot(SAL_int,Southinf_int,'r--')
	plot(SAL_int,Southsup_int,'r--')
	plot(Northinf_int,TEMP_int,'k--')
	plot(Northsup_int,TEMP_int,'k--')
pause;
end
%%%%%%%%% End Debug lines ....


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 			Historical Data			 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
grid_critSAF=zeros(size(grid_long));

if (CritSAF~=0)
for i=1:length(grid_long)
	isok=find(~isnan(grid_pres(:,i))&~isnan(grid_sal(:,i))&~isnan(grid_ptmp(:,i)) & [diff(grid_pres(:,i)); 0]~=0);
        if( isempty(find(diff(grid_pres(isok,i))==0)) & length(isok)>2 & min(grid_pres(isok,i))<300 & max(grid_pres(isok,i))>300)
		temp=sw_temp(grid_sal(isok,i), grid_ptmp(isok,i), grid_pres(isok,i), 0);
   		T300=interp1(grid_pres(isok,i),temp,300);
		if (T300>5) 
			grid_critSAF(i)=1;
		elseif (T300<3)
			grid_critSAF(i)=-1;
		else
			grid_critSAF(i)=0;
		end
	end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 			Test				 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
isSelect=find(grid_critSAF==CritSAF);

hist_sal	= grid_sal(:,isSelect);
hist_ptmp	= grid_ptmp(:,isSelect);
hist_pres	= grid_pres(:,isSelect);
hist_lat	= grid_lat(isSelect);
hist_lon	= grid_long(isSelect);
hist_dates	= grid_dates(isSelect);
hist_Z		= grid_Z(isSelect);


else %if Lat<-30
hist_sal	= grid_sal;
hist_ptmp	= grid_ptmp;
hist_pres	= grid_pres;
hist_lat	= grid_lat;
hist_lon	= grid_long;
hist_dates	= grid_dates;
hist_Z		= grid_Z;
end%if Lat<-30


return 
