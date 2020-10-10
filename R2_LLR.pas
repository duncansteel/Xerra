program R2_LLR(input, output, stdErr);

{Duncan Steel, 2020 October 06}

{Program for reading in Radarsat-2 ephemeris written out of STK, this time in Fixed LLR coordinates}
{Note latitude is GEOCENTRIC not GEODETIC}   {In STK, LLA implies geodetic}

var
  i,j,k,n,dd,mm,yy,hh,mn,ss : integer;
  r,z,p,cc,beta0,dbeta,beta,fbeta,fdbeta,omega,twobet,betmom,aob,h : real;
  sma,smb,oneof,gmearth,e2,fss,clat,dlat,lon,rad,phi,phid,pi,degrad,raddeg,sma2,smb2,dsm2 : real;
  radphi,sphi,cphi,dummy1,dummy2,dummy3,dummy4,dummya,dummyb,dummyc,dummyd : real;
  ifl,jfl,ofl,pfl,rfl : text;
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
aob:=sma/smb;
oneof:=298.257223563; {Reciprocal of WGS84 ellipsoid flattening}
e2:=1.0 - sqr(smb/sma); {Square of eccentricity of WGS84 ellipsoid}
dummy1:=1.0/(1.0 - e2);
dummy2:=1.0/oneof;
sma2:=sqr(sma);
smb2:=sqr(smb);
dsm2:=sma2 - smb2;

assign(ifl,'Radarsat-2_Fixed_LLR_Position.txt');
reset(ifl);
assign(ofl,'R2_LLR.txt');
rewrite(ofl);
assign(pfl,'R2_LLR.csv');
rewrite(pfl);
writeln; writeln;

for i:=1 to 6 do readln(ifl);
j:=-1; {counter}
k:=0; {another counter}
while not(eof(ifl)) do begin
read(ifl,dd);
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
read(ifl,clat); {geoCENTRIC latitude,deg}
read(ifl,lon); {longitude, deg}
read(ifl,rad); {radial distance from Earth centre, km}
readln(ifl);
if (abs(fss - ss) < 0.000000001) then begin {a change in orbital element set causes a part-second jump, so ignore it}
j:=j+1;
k:=k+1;
phid:=clat*degrad; {geocentric latitude in radians}
{Method to derive geodetic latitude dlat from Burtch (2006), Borkowski's interative method}
r:=rad*1000.0;
p:=r*cos(phid);
z:=r*sin(phid);
dummyc:=sma*p;
dummyd:=smb*z;
cc:=(dsm2)/(sqrt(sqr(dummyc) + sqr(dummyd)));
omega:=arctan(dummyd/dummyc);
beta0:=arctan((sma*z)/(smb*p));
beta:=beta0;
dbeta:=0.001;
while (abs(dbeta) > 0.000000001) do begin
twobet:=2.0*beta;
betmom:=(beta - omega);
fbeta:=(2.0*sin(betmom)) - (cc*sin(twobet));
fdbeta:=2.0*((cos(betmom)) - (cc*cos(twobet)));
dbeta:=fbeta/fdbeta;
beta:=beta - dbeta;
                             end; {while abs dbeta}
phi:=(arctan(aob*sin(beta)/cos(beta)));    {geodetic latitude in radians}
dlat:=phi*raddeg;           {in degrees}
h:=(((p - (sma*cos(beta)))*cos(phi)) + ((z - (smb*sin(beta)))*sin(phi)))/1000.0;   {geodetic altitude in km}
sphi:=sin(phi);
cphi:=cos(phi);
dummya:=sqr(sma2*cphi);
dummyb:=sqr(smb2*sphi);
dummyc:=sqr(sma*cphi);
dummyd:=sqr(smb*sphi);
radphi:=sqrt((dummya + dummyb) / (dummyc + dummyd))/1000.0; {Radius of Earth at this geodetic latitude, in km}
writeln(ofl,j+1:8,yy:6,mm:4,dd:4,hh:4,mn:4,ss:4,clat:10:3,lon:10:3,rad:14:6,radphi:14:6,h:14:6,dlat:12:6);
writeln(pfl,j+1:8,' , ',yy:6,' , ',mm:4,' , ',dd:4,' , ',hh:4,' , ',mn:4,' , ',ss:4,' , ',dlat:10:3,' , ',
            clat:10:3,' , ',lon:10:3,' , ',rad:12:6,' , ',radphi:12:6,' , ',h:12:6);
                                        end;    {abs(fss - ss) < 0.000000001}
end; {while not eof ifl}
close(ifl);
close(ofl);
close(pfl);
{Now test the latitudes}
assign(rfl,'R2_LLR.txt');
reset(rfl);
assign(jfl,'Radarsat-2_LLA.txt');
reset(jfl);
for i:=1 to k do begin  {Counter from above}
for j:=1 to 25 do read(jfl,c);
readln(jfl,dummya,dummyb,dummyc);  {geodetic latitude, longitude, and altitude}  {Ignore Cartesian coordinates}
for j:=1 to 46 do read(rfl,c);     {read to past the geocentric latitude}
read(rfl,dummy1);                  {longitude}
if (abs(dummy1 - dummyb) > 0.001) then writeln('Line',i:10,' - longitudes disagree');  {Compare longitudes}
read(rfl,dummy2,dummy3);  {radius vector, Earth radius at this latitude - not used}
read(rfl,dummy3); {geodetic height, my calculation}
readln(rfl,dummy4); {geodetic latitude, my calculation}
if (abs(dummy4 - dummya) > 0.001) then writeln('Line',i:10,' - latitudes disagree'); ;  {Compare geodetic latitudes}
if (abs(dummy3 - dummyc) > 0.001) then writeln('Line',i:10,' - heights disagree'); ;  {Compare geodetic heights}
if (i mod 1000) = 0 then writeln(i:10);
                   end; {i = 1 to k}
close(jfl); close(rfl);
readln;

end.
