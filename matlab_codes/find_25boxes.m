
% function [ pa_wmo_numbers ] = find_25boxes( pn_float_long, pn_float_lat, pa_wmo_boxes);
%
% This function finds 5x5=25 WMO boxes with the float profile in the centre.
% The WMO box numbers, between 90N and 90S, are stored in wmo_boxes.mat.
% The 1st column has the box numbers, the 2nd column denotes CTD data,
% the 3rd column denotes bottle data, the 4th column denotes Argo data.
% No data is denoted by 0. Otherwise 1.
%
% A. Wong, 16 August 2004
%
% C. Cabanes Nov. 2014 : extend la_x, so interp2 does not think longitudes in the range [5W 5E] are out-of-bound with matlab version >= R2012b

function [ pa_wmo_numbers ] = find_25boxes( pn_float_long, pn_float_lat, pa_wmo_boxes ) ;

pa_wmo_numbers = [ NaN.*ones( 25, 1 ), zeros( 25, 1 ), zeros( 25, 1 ), zeros( 25, 1 ) ] ;

la_lookup_x = [ ] ;
la_lookup_y = [ ] ;
vector_x = [] ;
vector_y = [] ;
%keyboard
la_x = [ -5:10:365 ] ; % 38 elements
for i=1:18
  la_lookup_x = [ la_lookup_x; la_x ] ;
end

la_y = [ 85:-10:-85 ] ; % 18 elements
for i=1:38
  la_lookup_y = [ la_lookup_y, la_y' ];
  vector_y = [ vector_y; la_y' ];
  vector_x = [ vector_x; la_x(i).*ones(18,1) ];
end

la_lookup_no=reshape( [ 1:648 ], 18, 36 ) ;
la_lookup_no=[la_lookup_no(:,end),la_lookup_no,la_lookup_no(:,1)];


ln_x1 = pn_float_long +.01;
ln_x2 = pn_float_long + 10.01;
ln_x3 = pn_float_long - 9.99;
ln_x4 = pn_float_long + 20.01;
ln_x5 = pn_float_long - 19.99;

ln_y1 = pn_float_lat + .01;
ln_y2 = pn_float_lat + 10.01;
ln_y3 = pn_float_lat - 9.99;
ln_y4 = pn_float_lat + 20.01;
ln_y5 = pn_float_lat - 19.99;

% interp2 will treat 360 as out of range, but will interpolate 0

if( ln_x3<0 )ln_x3=360+ln_x3;end;
if( ln_x5<0 )ln_x5=360+ln_x5;end;

if( ln_x1>=360 )ln_x1=ln_x1-360;end
if( ln_x2>=360 )ln_x2=ln_x2-360;end
if( ln_x4>=360 )ln_x4=ln_x4-360;end

if( isnan(pn_float_lat)==0&isnan(pn_float_long)==0 )
	ln_i1 = interp2( la_lookup_x, la_lookup_y, la_lookup_no, ln_x1, ln_y1, 'nearest' ) ;
	ln_i2 = interp2( la_lookup_x, la_lookup_y, la_lookup_no, ln_x2, ln_y1, 'nearest' ) ;
	ln_i3 = interp2( la_lookup_x, la_lookup_y, la_lookup_no, ln_x3, ln_y1, 'nearest' ) ;
	ln_i4 = interp2( la_lookup_x, la_lookup_y, la_lookup_no, ln_x4, ln_y1, 'nearest' ) ;
	ln_i5 = interp2( la_lookup_x, la_lookup_y, la_lookup_no, ln_x5, ln_y1, 'nearest' ) ;
	ln_i6 = interp2( la_lookup_x, la_lookup_y, la_lookup_no, ln_x1, ln_y2, 'nearest' ) ;
	ln_i7 = interp2( la_lookup_x, la_lookup_y, la_lookup_no, ln_x2, ln_y2, 'nearest' ) ;
	ln_i8 = interp2( la_lookup_x, la_lookup_y, la_lookup_no, ln_x3, ln_y2, 'nearest' ) ;
	ln_i9 = interp2( la_lookup_x, la_lookup_y, la_lookup_no, ln_x4, ln_y2, 'nearest' ) ;
	ln_i10= interp2( la_lookup_x, la_lookup_y, la_lookup_no, ln_x5, ln_y2, 'nearest' ) ;
	ln_i11= interp2( la_lookup_x, la_lookup_y, la_lookup_no, ln_x1, ln_y3, 'nearest' ) ;
	ln_i12= interp2( la_lookup_x, la_lookup_y, la_lookup_no, ln_x2, ln_y3, 'nearest' ) ;
	ln_i13= interp2( la_lookup_x, la_lookup_y, la_lookup_no, ln_x3, ln_y3, 'nearest' ) ;
	ln_i14= interp2( la_lookup_x, la_lookup_y, la_lookup_no, ln_x4, ln_y3, 'nearest' ) ;
	ln_i15= interp2( la_lookup_x, la_lookup_y, la_lookup_no, ln_x5, ln_y3, 'nearest' ) ;
	ln_i16= interp2( la_lookup_x, la_lookup_y, la_lookup_no, ln_x1, ln_y4, 'nearest' ) ;
	ln_i17= interp2( la_lookup_x, la_lookup_y, la_lookup_no, ln_x2, ln_y4, 'nearest' ) ;
	ln_i18= interp2( la_lookup_x, la_lookup_y, la_lookup_no, ln_x3, ln_y4, 'nearest' ) ;
	ln_i19= interp2( la_lookup_x, la_lookup_y, la_lookup_no, ln_x4, ln_y4, 'nearest' ) ;
	ln_i20= interp2( la_lookup_x, la_lookup_y, la_lookup_no, ln_x5, ln_y4, 'nearest' ) ;
	ln_i21= interp2( la_lookup_x, la_lookup_y, la_lookup_no, ln_x1, ln_y5, 'nearest' ) ;
	ln_i22= interp2( la_lookup_x, la_lookup_y, la_lookup_no, ln_x2, ln_y5, 'nearest' ) ;
	ln_i23= interp2( la_lookup_x, la_lookup_y, la_lookup_no, ln_x3, ln_y5, 'nearest' ) ;
	ln_i24= interp2( la_lookup_x, la_lookup_y, la_lookup_no, ln_x4, ln_y5, 'nearest' ) ;
	ln_i25= interp2( la_lookup_x, la_lookup_y, la_lookup_no, ln_x5, ln_y5, 'nearest' ) ;
else
        ln_i1 = NaN ;
        ln_i2 = NaN ;
        ln_i3 = NaN ;
        ln_i4 = NaN ;
        ln_i5 = NaN ;
        ln_i6 = NaN ;
        ln_i7 = NaN ;
        ln_i8 = NaN ;
        ln_i9 = NaN ;
        ln_i10 = NaN ;
        ln_i11 = NaN ;
        ln_i12 = NaN ;
        ln_i13 = NaN ;
        ln_i14 = NaN ;
        ln_i15 = NaN ;
        ln_i16 = NaN ;
        ln_i17 = NaN ;
        ln_i18 = NaN ;
        ln_i19 = NaN ;
        ln_i20 = NaN ;
        ln_i21 = NaN ;
        ln_i22 = NaN ;
        ln_i23 = NaN ;
        ln_i24 = NaN ;
        ln_i25 = NaN ;
end


if( isnan(ln_i1)==0 )pa_wmo_numbers(1,:)=pa_wmo_boxes(ln_i1,:);end;
if( isnan(ln_i2)==0 )pa_wmo_numbers(2,:)=pa_wmo_boxes(ln_i2,:);end;
if( isnan(ln_i3)==0 )pa_wmo_numbers(3,:)=pa_wmo_boxes(ln_i3,:);end;
if( isnan(ln_i4)==0 )pa_wmo_numbers(4,:)=pa_wmo_boxes(ln_i4,:);end;
if( isnan(ln_i5)==0 )pa_wmo_numbers(5,:)=pa_wmo_boxes(ln_i5,:);end;
if( isnan(ln_i6)==0 )pa_wmo_numbers(6,:)=pa_wmo_boxes(ln_i6,:);end;
if( isnan(ln_i7)==0 )pa_wmo_numbers(7,:)=pa_wmo_boxes(ln_i7,:);end;
if( isnan(ln_i8)==0 )pa_wmo_numbers(8,:)=pa_wmo_boxes(ln_i8,:);end;
if( isnan(ln_i9)==0 )pa_wmo_numbers(9,:)=pa_wmo_boxes(ln_i9,:);end;
if( isnan(ln_i10)==0 )pa_wmo_numbers(10,:)=pa_wmo_boxes(ln_i10,:);end;
if( isnan(ln_i11)==0 )pa_wmo_numbers(11,:)=pa_wmo_boxes(ln_i11,:);end;
if( isnan(ln_i12)==0 )pa_wmo_numbers(12,:)=pa_wmo_boxes(ln_i12,:);end;
if( isnan(ln_i13)==0 )pa_wmo_numbers(13,:)=pa_wmo_boxes(ln_i13,:);end;
if( isnan(ln_i14)==0 )pa_wmo_numbers(14,:)=pa_wmo_boxes(ln_i14,:);end;
if( isnan(ln_i15)==0 )pa_wmo_numbers(15,:)=pa_wmo_boxes(ln_i15,:);end;
if( isnan(ln_i16)==0 )pa_wmo_numbers(16,:)=pa_wmo_boxes(ln_i16,:);end;
if( isnan(ln_i17)==0 )pa_wmo_numbers(17,:)=pa_wmo_boxes(ln_i17,:);end;
if( isnan(ln_i18)==0 )pa_wmo_numbers(18,:)=pa_wmo_boxes(ln_i18,:);end;
if( isnan(ln_i19)==0 )pa_wmo_numbers(19,:)=pa_wmo_boxes(ln_i19,:);end;
if( isnan(ln_i20)==0 )pa_wmo_numbers(20,:)=pa_wmo_boxes(ln_i20,:);end;
if( isnan(ln_i21)==0 )pa_wmo_numbers(21,:)=pa_wmo_boxes(ln_i21,:);end;
if( isnan(ln_i22)==0 )pa_wmo_numbers(22,:)=pa_wmo_boxes(ln_i22,:);end;
if( isnan(ln_i23)==0 )pa_wmo_numbers(23,:)=pa_wmo_boxes(ln_i23,:);end;
if( isnan(ln_i24)==0 )pa_wmo_numbers(24,:)=pa_wmo_boxes(ln_i24,:);end;
if( isnan(ln_i25)==0 )pa_wmo_numbers(25,:)=pa_wmo_boxes(ln_i25,:);end;
%keyboard
