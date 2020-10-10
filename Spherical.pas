program Spherical(input, output);

{Duncan Steel, 2020 September 28}
{Output folder changed 2020 October 06}

var
  i : integer;
  pi,piby2,twopi,degrad,raddeg,sma,smb,sma2,smb2,la1,la2,lo1,lo2,dlat,dlon,theta,d,t1,t2,t3,t4,r,psi,
       sdlon,cdlon,sp1,cp1,sp2,cp2,st,ct,alp1,alp2,deltsig : real;
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

function atan2(sa,ca : real) : real;
{Returns the inverse tangent of an angle a where sa=sin(a) and ca=cos(a) }
var chk : real;
begin
chk:=1.0e30;
if ((abs(sa)<1.0e-12) and (abs(ca)<1.0e-12)) then begin      {both parameters are zero}
          writeln('Error 1 in atan2 function'); chk:=0.0;
                                                    end;
if ((abs(ca)<1.0e-12) and (sa>0.0)) then chk:=piby2;
if ((abs(ca)<1.0e-12) and (sa<0.0)) then chk:=-piby2;
if (ca>0.0) then chk:=arctan(sa/ca);
if ((ca<0.0) and (sa>=0.0)) then chk:=(arctan(sa/ca)) + pi;
if ((ca<0.0) and (sa<0.0)) then chk:=(arctan(sa/ca)) - pi;
if (chk>9.9e29) then writeln('Error 2 in atan2 function');
if (chk < 0.0) then chk:=chk + twopi;
atan2:=chk;
end; {of function atan2}

{------------------------------------------------------------------}

begin
assign(ofl,'Spherical.txt');
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
la1:=41.0;  {Point P1}
lo2:=40.0;
la2:=40.0;  {Point P2}
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
dlat:=la2-la1;
dlon:=lo2-lo1;
sp1:=sin(la1);
cp1:=cos(la1);
sp2:=sin(la2);
cp2:=cos(la2);
sdlon:=sin(dlon);
cdlon:=cos(dlon);
{writeln(sdlon:10:5,cdlon:10:5,sp1:10:5,cp1:10:5,sp2:10:5,cp2:10:5); }
psi:=(la1+la2)/2.0;
t1:=sin(psi);
t2:=cos(psi);
t3:=sqr(sma2*t2) + sqr(smb2*t1);
t4:=sqr(sma*t2) + sqr(smb*t1);
r:=sqrt(t3/t4); {Earth radius at average latitude}

t1:=cp2*sdlon;
t2:=(cp1*sp2) - (sp1*cp2*cdlon);
alp1:=atan2(t1,t2);

st:=sqrt(sqr(t2) + sqr(t1));
ct:=(sp1*sp2) + (cp1*cp2*cdlon);
deltsig:=atan2(st,ct);

t1:=cp1*sdlon;
t2:=(sp2*cp1*cdlon) - (cp2*sp1);
alp2:=atan2(t1,t2);

d:=r*deltsig;

writeln(ofl);
writeln(ofl,alp1*raddeg:12:6,alp2*raddeg:12:6,deltsig*raddeg:12:6,d/1000.0:15:3);


writeln(ofl);
writeln('Done - hit return'); writeln;

readln;
close(ofl);

end.
