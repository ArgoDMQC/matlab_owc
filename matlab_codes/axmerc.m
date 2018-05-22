
ax = axis;
axis square;


delx = ax(2) - ax(1);
dely = ax(4) - ax(3);
latmid = (ax(3)+ax(4))/2;

scaley = delx*cos(pi*latmid/180)/dely;

if ( scaley > 1)
% keep lon, scale lat

ax(3) = latmid - delx*cos(pi*latmid/180)/2;
ax(4) = latmid + delx*cos(pi*latmid/180)/2;

end

if ( scaley < 1)
%keep lat scale lon

lonmid = ( ax(1) + ax(2) ) /2;

ax(1) = lonmid - dely/cos(pi*latmid/180)/2;
ax(2) = lonmid + dely/cos(pi*latmid/180)/2;

end


axis(ax);
