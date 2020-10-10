program Haversinsteps(input, output);

{Duncan Steel, 2020 September 28}
{Output folder changed 2020 October 06}

{To find a path from point A to point B defined by lat & lon, as a great circle defined by haversines
 and 1000 steps in longitude, each section being a rhumb line}
{This is for comparison with Vincenty calculations}

var
  i : integer;
  pi,piby2,twopi,degrad,raddeg,sma,smb,sma2,smb2,la1,la2,lo1,lo2,dlat,dlon,theta,d,t1,t2,t3,t4,r,psi : real;
  ofl : text;
  c : char;


{------------------------------------------------------------------}

FUNCTION ASIN(XYZ : REAL) : REAL;

{Returns the inverse sine of parameter XYZ}

BEGIN
IF ABS(XYZ) > 0.999999 THEN ASIN:=XYZ*PI/2.0
     ELSE ASIN:=ARCTAN(XYZ/SQRT(1.0-(XYZ*XYZ)));
END;   {of function ASIN}

{------------------------------------------------------------------}

begin
assign(ofl,'Haversinsteps.txt');
rewrite(ofl);
writeln; writeln;

c:=' ';
pi:=4.0*arctan(1.0);
piby2:=pi/2.0;
twopi:=2.0*pi;
degrad:=pi/180.0;
raddeg:=180.0/pi;
sma:=6378137.0;    {Semi-major axis of WGS84 ellipsoid}
smb:=6356752.314245; {Semi-minor axis of WGS84 ellipsoid}
sma2:=sqr(sma);
smb2:=sqr(smb);

lo1:=41.0; {These are the overall start and end points: longitude and latitude}
la1:=41.0;  {Point A}
lo2:=40.0;
la2:=40.0;  {Point B}

writeln(ofl);

if abs(la1) > 90.0 then writeln('Error in initial latitude');
if abs(la2) > 90.0 then writeln('Error in final latitude');
if abs(lo1) > 360.0 then writeln('Error in initial longitude');
if abs(lo2) > 360.0 then writeln('Error in final longitude');
while (lo1 < 0.0) do lo1:=lo1 + 360.0;
while (lo2 < 0.0) do lo2:=lo2 + 360.0;
lo1:=lo1*degrad;
la1:=la1*degrad;
lo2:=lo2*degrad;
la2:=la2*degrad;
dlat:=la1-la2;
dlon:=lo1-lo2;
t1:=sqr(sin(dlat/2.0));
t2:=cos(la1)*cos(la2);
t3:=sqr(sin(dlon/2.0));
t4:=sqrt(t1 + (t2*t3));
theta:=2.0 * asin(t4);
psi:=(la1+la2)/2.0;
t1:=sin(psi);
t2:=cos(psi);
t3:=sqr(sma2*t2) + sqr(smb2*t1);
t4:=sqr(sma*t2) + sqr(smb*t1);
r:=sqrt(t3/t4);
d:=r * theta;
writeln(theta:12:6,theta*raddeg:12:6,r/1000.0:12:3,d/1000.0:15:6);
writeln(ofl,theta:12:6,theta*raddeg:12:6,r/1000.0:12:3,d/1000.0:15:6);

writeln(ofl);
writeln('Done - hit return'); writeln;

readln;
close(ofl);

end.
