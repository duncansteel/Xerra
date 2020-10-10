program WayPoints(input, output);

{Duncan Steel, 2020 October 06}
{Output folder changed 2020 October 06}

var
  i,j,n : integer;
  pi,piby2,twopi,degrad,raddeg,sma,smb,sma2,smb2,la1,la2,lo1,lo2,dlat,dlon,theta,d,t1,t2,t3,t4,r,psi,
       sdlon,cdlon,sp1,cp1,sp2,cp2,st,ct,alp1,alp2,alp0,deltsig,sig01,sig02,lam0,lam01,cotalp0,talp0,
       salp0,calp0,loni,lati,sigi,llon,llat,dlonn,lsig,totd,xd,alp,x,y : real;
  ofl,cfl,afl,rfl : text;
  c : char;

{------------------------------------------------------------------}

FUNCTION ASIN(XYZ : REAL) : REAL;

{Returns the inverse sine of parameter XYZ}

BEGIN
IF ABS(XYZ) > 0.999999 THEN ASIN:=XYZ*PI/2.0
     ELSE ASIN:=ARCTAN(XYZ/SQRT(1.0-(XYZ*XYZ)));
END;   {of function ASIN}

{------------------------------------------------------------------}

FUNCTION ACOS(ABC : REAL) : REAL;

{Returns the inverse cosine of parameter ABC}

BEGIN
IF ABS(ABC) > 0.999999 THEN BEGIN
                            IF ABC<0.0 THEN ACOS:=PI
                                       ELSE ACOS:=0.0;
                              END
    ELSE ACOS:=(PI/2.0) - ARCTAN(ABC/SQRT(1.0-(ABC*ABC)));
END; {of function ACOS}

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
assign(ofl,'WayPoints.txt');
rewrite(ofl);
assign(cfl,'WayPoints.csv');
rewrite(cfl);
assign(afl,'WayPointsAz.csv');
rewrite(afl);
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

lo1:=30.0; {These are the overall start and end points: longitude and latitude}
la1:=35.0; {Point X1}
lo2:=29.0;
la2:=34.0; {Point X2}

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

t1:=cp2*sdlon;
t2:=(cp1*sp2) - (sp1*cp2*cdlon);
alp1:=atan2(t1,t2);

t3:=cp1*sdlon;
t4:=(sp2*cp1*cdlon) - (cp2*sp1);
alp2:=atan2(t3,t4);

t3:=sin(alp1);
t4:=cos(alp1);
st:=t3*cp1;
ct:=sqrt(sqr(t4) + sqr(t3*sp1));
alp0:=atan2(st,ct);
talp0:=st/ct;
salp0:=st;
calp0:=ct;
cotalp0:=1.0/talp0;

st:=sp1/cp1;
ct:=t4;
sig01:=atan2(st,ct);

st:=sp2/cp2;
ct:=cos(alp2);
sig02:=atan2(st,ct);

st:=sqrt(sqr(t2) + sqr(t1));
ct:=(sp1*sp2) + (cp1*cp2*cdlon);
deltsig:=atan2(st,ct);

psi:=(la1+la2)/2.0; {average latitude, to get radius r}
t1:=sin(psi);
t2:=cos(psi);
t3:=sqr(sma2*t2) + sqr(smb2*t1);
t4:=sqr(sma*t2) + sqr(smb*t1);
r:=sqrt(t3/t4); {Earth radius at average latitude}

{
n:=round(abs(psi*raddeg*100.0));
assign(rfl,'EarthRadii.txt');
reset(rfl);
for i:=0 to n do readln(rfl,j,y,x);
close(rfl);
}

xd:=r*deltsig; {total distance moved in metres from point 1 to point 2}
writeln(ofl);
writeln(ofl,'Start azimuth  End azimuth   Distance (km) ');
writeln(ofl,alp1*raddeg:12:6,alp2*raddeg:12:6,xd/1000.0:15:3);
writeln(ofl);
writeln(ofl,'  Step   Longitude   Latitude        Azimuth      Step length (km)');
i:=1;
writeln(afl,i:6,' , ',alp1*raddeg:10:6);

st:=sin(alp0) * sin(sig01);
ct:=cos(sig01);
lam01:=atan2(st,ct);
lam0:=lo1 - lam01;
if (lam0 < 0.0) then lam0:=lam0 + twopi;

n:=100;
totd:=0.0;
llon:=lo1;
llat:=la1;
lsig:=sig01;
dlonn:=dlon/n;
for i:=1 to n do begin
loni:=llon + dlonn;
t2:=loni - lam0;
t3:=sin(t2);
t4:=cos(t2);
lati:=arctan(cotalp0 * t3);
sigi:=acos(cos(lati) * t4);
deltsig:=sigi - lsig;
t1:=cos(sigi);
alp:=atan2(talp0,t1);
writeln(afl,i+1:6,' , ',alp*raddeg:10:6);

psi:=(llat+lati)/2.0; {average latitude, to get radius r}
t1:=sin(psi);
t2:=cos(psi);
t3:=sqr(sma2*t2) + sqr(smb2*t1);
t4:=sqr(sma*t2) + sqr(smb*t1);
r:=sqrt(t3/t4); {Earth radius at average latitude}
d:=r*deltsig; {distance moved in metres}
totd:=totd + d;
writeln(ofl,i:6,loni*raddeg:12:6,lati*raddeg:12:6,alp*raddeg:16:6,'  ',d/1000.0:12:6);
writeln(cfl,i:6,' , ',llat*raddeg:10:6,' , ',llon*raddeg:10:6);
llat:=lati;
llon:=loni;
lsig:=sigi;
                     end; {i=1 to n}
writeln(ofl,'                            Total distance = ',totd/1000.0:12:3,' km');
writeln(cfl,n+1:6,' , ',llat*raddeg:10:6,' , ',llon*raddeg:10:6);

writeln(ofl);
close(ofl);
close(cfl);
close(afl);

end.
