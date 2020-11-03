program AT3(input, output, stdErr);

{Duncan Steel, 2020 October 28}

{AT = Auto Tasking}

{Based on R2_LLR.pas...}
{Program for reading in Radarsat-2 ephemeris written out of STK, this time in Fixed LLR coordinates}
{Note latitude is GEOCENTRIC not GEODETIC}   {In STK, LLA implies geodetic}

var
  j,k,l : LongInt;
  i,n,day,mm,yy,hh,mn,ss : integer;
  r,ctarlat,tarlat,tarlon,xt,yt,zt,xs,ys,zs,slat,clat,sma,smb,gmearth,e2,fss,lon,rlim,happ,minr,
  satrad,pi,degrad,raddeg,sma2,smb2,tarrad,satlat,satlon,b2oa2,dum,t1,t2,t3,t4,dx,dy,dz : real;
  ifl,ofl,pfl : text;
  m : array[1..3] of char;
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

begin
c:=' ';
for i:=1 to 3 do m[i]:=c;
pi:=4.0*arctan(1.0);
degrad:=pi/180.0;
raddeg:=180.0/pi;
gmearth:=3.9860044188E14; {gravitational constant * mass of the Earth in SI units}
sma:=6378137.0;    {Semi-major axis of WGS84 ellipsoid}
smb:=6356752.314245; {Semi-minor axis of WGS84 ellipsoid}
sma2:=sqr(sma);
smb2:=sqr(smb);
e2:=1.0 - (smb2/sma2); {Square of eccentricity of WGS84 ellipsoid}
b2oa2:=smb2/sma2;

writeln; writeln;
{
writeln('Target geodetic latitude? '); readln(tarlat);
writeln;
writeln('Target longitude?         '); readln(tarlon);
writeln;
}
tarlat:=-23.6417;
tarlon:=-178.911; {North Minerva Reef}

tarlat:=tarlat*degrad;
tarlon:=tarlon*degrad;
{On the surface so h or A = 0}

ctarlat:=arctan(b2oa2*sin(tarlat)/cos(tarlat)); {geocentric latitude of the target in radians}
slat:=sin(ctarlat);
clat:=cos(ctarlat);

t1:=sin(tarlat);
t2:=cos(tarlat);
t3:=sqr(sma2*t2) + sqr(smb2*t1);
t4:=sqr(sma*t2) + sqr(smb*t1);
tarrad:=sqrt(t3/t4); {Earth radius at geodetic latitude tarlat}

xt:=tarrad*clat*cos(tarlon);
yt:=tarrad*clat*sin(tarlon);
zt:=tarrad * slat;

assign(ifl,'R2_ephem.txt');
reset(ifl);
assign(ofl,'AT3LLRdelta-r.txt');
rewrite(ofl);
assign(pfl,'AT1xyx.txt');
rewrite(pfl);
writeln; writeln;
writeln(ofl);
writeln(ofl,'Target geocentric latitude, longitude, radius: ',
              ctarlat*raddeg:12:4,tarlon*raddeg:12:4,tarrad/1000.0:12:4);
tarrad:=tarrad/1000.0; {in km}
writeln(ofl);

for i:=1 to 6 do readln(ifl);
j:=-1; {counter}
k:=0;
minr:=2500.0; {initial value for minimum distance/range}

while not(eof(ifl)) do begin
read(ifl,day);
read(ifl,c);
for i:=1 to 3 do begin
read(ifl,c);
m[i]:=c;
                   end;
mm:=0;
if ((m[1] = 'J') and (m[2] = 'a') and (m[3] = 'n')) then mm:=1;
if ((m[1] = 'F') and (m[2] = 'e') and (m[3] = 'b')) then mm:=2;
if ((m[1] = 'M') and (m[2] = 'a') and (m[3] = 'r')) then mm:=3;
if ((m[1] = 'A') and (m[2] = 'p') and (m[3] = 'r')) then mm:=4;
if ((m[1] = 'M') and (m[2] = 'a') and (m[3] = 'y')) then mm:=5;
if ((m[1] = 'J') and (m[2] = 'u') and (m[3] = 'n')) then mm:=6;
if ((m[1] = 'J') and (m[2] = 'u') and (m[3] = 'l')) then mm:=7;
if ((m[1] = 'A') and (m[2] = 'u') and (m[3] = 'g')) then mm:=8;
if ((m[1] = 'S') and (m[2] = 'e') and (m[3] = 'p')) then mm:=9;
if ((m[1] = 'O') and (m[2] = 'c') and (m[3] = 't')) then mm:=10;
if ((m[1] = 'N') and (m[2] = 'o') and (m[3] = 'v')) then mm:=11;
if ((m[1] = 'D') and (m[2] = 'e') and (m[3] = 'c')) then mm:=12;
if (mm = 0) then writeln('Month lookup error');

read(ifl,yy);
read(ifl,c); {blank}
hh:=0;
read(ifl,c);
hh:=(ord(c) - 48) * 10;
read(ifl,c);
hh:=hh + (ord(c) - 48);
read(ifl,c);
mn:=0;
read(ifl,c);
mn:=(ord(c) - 48) * 10;
read(ifl,c);
mn:=mn + (ord(c) - 48);
read(ifl,c);
read(ifl,fss);
ss:=round(fss);
read(ifl,satlat); {satellite geoCENTRIC latitude,deg}
read(ifl,satlon); {satellite longitude, deg}
read(ifl,satrad); {radial distance from Earth centre, km}
readln(ifl);
if (abs(fss - ss) < 0.000000001) then begin {a change in orbital element set causes a part-second jump, so ignore it}
j:=j+1; {counter}
{if ((j mod 100)= 0) then writeln(j:6);}
satlat:=satlat*degrad; {geocentric latitude in radians}
satlon:=satlon*degrad; {longitude in radians}

dum:=cos(satlat);
xs:=satrad * dum * cos(satlon) * 1000.0;   {all in metres}
ys:=satrad * dum * sin(satlon) * 1000.0;
zs:=satrad * sin(satlat) * 1000.0;       {This gives answers in accord with STK to 0.1 m}

dx:=xs-xt;
dy:=ys-yt;
dz:=zs-zt;
r:=sqrt(sqr(dx) + sqr(dy) + sqr(dz)); {distance satellite to target, in metres}
  {This gives answers in accord with STK to within 2 or 3 metres}
r:=r/1000.0; {in km}
happ:=satrad - tarrad; {approximate height of satellite when near target}
rlim:=3.0 * happ; {a limit on the range, based on a slant angle from the satellite of 70 deg, and cos(70)=0.342}
if (r < rlim) then begin
if (r < minr) then begin
minr:=r;
k:=j;
                     end; {r < minr}
                     end; {r < rlim}
satlat:=satlat*raddeg;
satlon:=satlon*raddeg;
writeln(ofl,j:6,yy:6,mm:4,day:4,hh:4,mn:4,ss:4,satlat:14:6,satlon:14:6,satrad:14:6,r:14:6);
writeln(pfl,j:6,yy:6,mm:4,day:4,hh:4,mn:4,ss:4,xs/1000.0:15:6,ys/1000.0:15:6,zs/1000.0:15:6);
                                        end;    {abs(fss - ss) < 0.000000001}
end; {while not eof ifl}
close(ifl);
close(ofl);
close(pfl);

writeln(k:10,minr:12:3); readln;

assign(ifl,'AT3LLRdelta-r.txt');
reset(ifl);
assign(ofl,'MinR.txt');
rewrite(ofl);
for i:=1 to 3 do readln(ifl);
for j:=0 to (k-1) do readln(ifl,l);
readln(ifl,l,yy,mm,day,hh,mn,ss,satlat,satlon,satrad,r);
writeln(ofl,yy:6,mm:4,day:4,hh:4,mn:4,ss:4,satlat:14:6,satlon:14:6,satrad:14:6,r:14:6);
close(ifl);
close(ofl);

readln;

end.
