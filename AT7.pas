program AT7(input, output, stdErr);

{Duncan Steel, 2020 November 03}

{AT = Auto Tasking}

{Program for reading in a satellite ephemeris written out of STK in Fixed LLR coordinates,
 and then identifying overpasses near some terrestrial target}
{Note latitude is GEOCENTRIC not GEODETIC}   {In STK, LLA implies geodetic}

var
  j,k,dcnt : LongInt;
  i,n,day,mm,yy,hh,mn,ss,i1,i2,i3,i4,i5,i6,satno,loc,lday : integer;
  r,ctarlat,tarlat,tarlon,xt,yt,zt,xs,ys,zs,slat,clat,slon,clon,sma,smb,gmearth,e2,fss,lon,rlast,azim,elev,
  pi,piby2,twopi,degrad,raddeg,sma2,smb2,tarrad,satlat,satlon,satrad,b2oa2,zamin,zamax,satlatlast,
  dum,t1,t2,t3,t4,dx,dy,dz,dlat,dlon,rs,re,rz,ps,pe,pz,look,lazim,za,nadira,sslat,cslat,sslon,cslon : real;
  ifl,ofl,pfl,xfl : text;
  m : array[1..3] of char;
  s : array[1..50] of char;
  c,cc,sdir,shand : char;

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
c:=' ';
sdir:=' ';
for i:=1 to 3 do m[i]:=c;
for i:=1 to 50 do s[i]:=c;
pi:=4.0*arctan(1.0);
piby2:=pi/2.0;
twopi:=pi*2.0;
degrad:=pi/180.0;
raddeg:=180.0/pi;
gmearth:=3.9860044188E14; {gravitational constant * mass of the Earth in SI units}
sma:=6378137.0;    {Semi-major axis of WGS84 ellipsoid}
smb:=6356752.314245; {Semi-minor axis of WGS84 ellipsoid}
sma2:=sqr(sma);
smb2:=sqr(smb);
e2:=1.0 - (smb2/sma2); {Square of eccentricity of WGS84 ellipsoid}
b2oa2:=smb2/sma2;

satno:=0;
while ((satno < 1) or (satno > 7)) do begin
writeln; writeln;
writeln('Which satellite?   ');
writeln('1 = Radarsat-2');
writeln('2 = TerraSAR-X');
writeln('3 = TanDEM-X');
writeln('4 = Capella-Sequoia');
writeln('5 = NovaSAR-1');
writeln('6 = ICEYE-6');
writeln('7 = ICEYE-7');
write('Enter 1 - 7 :   ');
read(satno);
                                      end; {while}
writeln; writeln;

case satno of
1: assign(ifl,'Radarsat-2_ephem.txt');
2: assign(ifl,'TerraSAR-X_ephem.txt');
3: assign(ifl,'TanDEM-X_ephem.txt');
4: assign(ifl,'Capella-2_ephem.txt');
5: assign(ifl,'NovaSAR-1_ephem.txt');
6: assign(ifl,'ICEYE-X6_ephem.txt');
7: assign(ifl,'ICEYE-X7_ephem.txt');
end; {case}
reset(ifl);

loc:=-1;
while ((loc < 0) or (loc > 15)) do begin
writeln; writeln;
writeln('Which target location?'   );
writeln('0 = Some other location');
writeln('1 = North Minerva Reef');
writeln('2 = South Minerva Reef');
writeln('3 = Opua');
writeln('4 = Auckland Harbour');
writeln('5 = Tauranga');
writeln('6 = Wellington');
writeln('7 = Gisborne');
writeln('8 = New Plymouth');
writeln('9 = Nelson');
writeln('10 = Picton');
writeln('11 = Lyttelton');
writeln('12 = Dunedin');
writeln('13 = Bluff');
writeln('14 = Oban/Halfmoon Bay');
writeln('15 = Greymouth');
write('Enter 0 - 15:   ');
read(loc);
                                      end; {while}
writeln; writeln;

case loc of
0: begin
tarlat:=-100.0;
tarlon:=-200.0;
while ((abs(tarlat) > 90.0) or (abs(tarlon) > 180.0)) do begin
writeln('Enter geodetic/geographic coordinates as XX.XXXX; latitude between -90 and +90, longitude -180 to + 180 ');
writeln;
write('Enter latitude:  '); read(tarlat); writeln;
write('Enter longitude: '); read(tarlon); writeln;
                                                            end; {while}
     end; {Some other location}
1: begin tarlat:=-23.6417; tarlon:=-178.911; end; {North Minerva Reef}
2: begin tarlat:=-23.94; tarlon:=-179.13; end; {South Minerva Reef}
3: begin tarlat:=-35.31; tarlon:=174.12; end; {Opua}
4: begin tarlat:=-36.84; tarlon:=174.77; end; {Auckland}
5: begin tarlat:=-37.67; tarlon:=176.18; end; {Tauranga}
6: begin tarlat:=-41.29; tarlon:=174.78; end; {Wellington}
7: begin tarlat:=-38.67; tarlon:=178.03; end; {Gisborne}
8: begin tarlat:=-39.06; tarlon:=174.05; end; {New Plymouth}
9: begin tarlat:=-41.26; tarlon:=173.28; end; {Nelson}
10: begin tarlat:=-41.29; tarlon:=174.01; end; {Picton}
11: begin tarlat:=-43.61; tarlon:=172.71; end; {Lyttelton}
12: begin tarlat:=-45.87; tarlon:=170.53; end; {Dunedin}
13: begin tarlat:=-46.60; tarlon:=168.34; end; {Bluff}
14: begin tarlat:=-46.90; tarlon:=168.13; end; {Oban}
15: begin tarlat:=-42.45; tarlon:=171.20; end; {Greymouth}
end; {case}

tarlat:=tarlat*degrad;
tarlon:=tarlon*degrad;
{On the surface so h or A = 0}

zamin:=90.0;
zamax:=0.0;
while (zamin > zamax) do begin
i:=-1;
writeln; writeln;
writeln('What range of incident/zenith angles for beam at the target location?   ');
while ((i < 0) or (i > 80)) do begin
write('Enter minimum (degrees)   ');
read(i);
zamin:=i*1.0;
                                      end; {while}
writeln;
i:=-1;
while ((i < 1) or (i > 85)) do begin
write('Enter maximum (degrees)   ');
read(i);
zamax:=i*1.0;
                                      end; {while}
writeln; writeln;
if (zamin > zamax) then writeln('Error in input of min/max angles');
writeln; writeln;
end; {while zamin > zamax}

ctarlat:=arctan(b2oa2*sin(tarlat)/cos(tarlat)); {geocentric latitude of the target in radians}
slat:=sin(ctarlat);
clat:=cos(ctarlat);
slon:=sin(tarlon);
clon:=cos(tarlon);

t1:=sin(tarlat);
t2:=cos(tarlat);
t3:=sqr(sma2*t2) + sqr(smb2*t1);
t4:=sqr(sma*t2) + sqr(smb*t1);
tarrad:=sqrt(t3/t4); {Earth radius at geodetic latitude tarlat}

xt:=tarrad*clat*cos(tarlon);
yt:=tarrad*clat*sin(tarlon);
zt:=tarrad * slat;
tarrad:=tarrad/1000.0; {in km}

assign(ofl,'LLRdr.txt');
rewrite(ofl);
assign(pfl,'EphemXYZ.txt');
rewrite(pfl);

readln(ifl);
i:=1;
while ((c <> ':') and (i < 50)) do begin
   read(ifl,c);
   s[i]:=c;
   i:=i+1;
                      end;
readln(ifl);
for n:=1 to 4 do readln(ifl);

ctarlat:=ctarlat*raddeg;
tarlon:=tarlon*raddeg;
write(ofl,'Target geocentric latitude, longitude, radius:',ctarlat:9:4,tarlon:10:4,tarrad:11:4);
write(ofl,'  Target geodetic latitude:',tarlat*raddeg:9:4);
writeln(ofl,'  X,Y,Z: ',xt/1000.0:14:6,yt/1000.0:14:6,zt/1000.0:14:6);
for n:=1 to (i-2) do write(ofl,s[n]);
writeln(ofl);
write(ofl,'  Index    Year   Mon   Day  Hour   Min   Sec   ');
writeln(ofl,'Geocentric Lat         Lon         Distance     Target Distance');

j:=-1; {counter}
dcnt:=0;
lday:=0;
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
if (yy < 2021) then begin {only to the end of 2020}
dcnt:=dcnt+1;
if ((dcnt mod 8640) = 1) then begin   {This is for ephemerides with 10-second jumps}
writeln('Working on day number ',((dcnt div 8640)+1):4)
                                 end;

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
satlat:=satlat*raddeg;
satlon:=satlon*raddeg;
writeln(ofl,j:7,yy:8,mm:6,day:6,hh:6,mn:6,ss:6,satlat:16:6,satlon:16:6,satrad:16:6,r:16:6);
writeln(pfl,j:6,yy:6,mm:4,day:4,hh:4,mn:4,ss:4,xs/1000.0:15:6,ys/1000.0:15:6,zs/1000.0:15:6);
                                        end;    {abs(fss - ss) < 0.000000001}

end {if yy<2021} else begin
readln(ifl);
if (day <> lday) then begin write('Tick...'); if (day mod 5) = 0 then writeln; end;
lday:=day;
                        end; {else}
end; {while not eof ifl}
close(ifl);
close(ofl);
close(pfl);

assign(ifl,'LLRdr.txt');
reset(ifl);
assign(ofl,'MinimumDistances.txt');
rewrite(ofl);
for i:=1 to 3 do readln(ifl);

rlast:=3100.0;
n:=0;
j:=0;
satlatlast:=-90.0;
while not(eof(ifl)) do begin
readln(ifl,j,yy,mm,day,hh,mn,ss,satlat,satlon,satrad,r);
if (satlat > satlatlast) then sdir:='A' else sdir:='D'; {Satellite ascending or descending (node) }
satlatlast:=satlat;
if (j > 0) then begin
if (r < 3000.0) then begin
if (r < rlast) then begin
rlast:=r;
                      end {r < rlast}
                      else begin
                      if (n = 0) then
          begin writeln(ofl,(j-1):6,i1:6,i2:4,i3:4,i4:4,i5:4,i6:4,t1:14:6,t2:14:6,t3:14:6,t4:14:6,'  ',sdir);
                        n:=1;
                        end;
                             end; {else}
                       end {r < 3000} else begin rlast:=3100.0; n:=0; end; {else}
                  end; {j > 0}
i1:=yy;
i2:=mm;
i3:=day;
i4:=hh;
i5:=mn;
i6:=ss;
t1:=satlat;
t2:=satlon;
t3:=satrad;
t4:=r;
                         end; {while not eof}
close(ifl);
close(ofl);

assign(ifl,'LLRdr.txt');
reset(ifl);
c:='X';
while (c <> 'Z') do read(ifl,c);
read(ifl,c);
read(ifl,c);
read(ifl,xt,yt,zt);
close(ifl);
writeln;
j:=0;
k:=-1;
n:=0;
case satno of
1: assign(ofl,'Radarsat-2_opportunities.txt');
2: assign(ofl,'TerraSAR-X_opportunities.txt');
3: assign(ofl,'TanDEM-X_opportunities.txt');
4: assign(ofl,'Capella-2_opportunities.txt');
5: assign(ofl,'NovaSAR-1_opportunities.txt');
6: assign(ofl,'ICEYE-6_opportunities.txt');
7: assign(ofl,'ICEYE-7_opportunities.txt');
end; {case}
rewrite(ofl);

assign(xfl,'LLRdr.txt');
reset(xfl);
n:=0;
c:=' ';
while (c <> 'X') do begin
read(xfl,c);
n:=n+1;
if (c <> 'X') then begin
write(ofl,c);        end;
                            end;    {while}
writeln(ofl);
readln(xfl,c);
c:=' ';
for n:=1 to 50 do s[n]:=c;
n:=0;
while (not(eoln(xfl))) do begin
read(xfl,c);
n:=n+1;
if (n = 10) then c:=' ';
s[n]:=c;
write(ofl,c);
                            end;
case satno of
1: write(ofl,'   - Note MDA slew plan provides for beam Right only, through to end-May 2021 ');
2: write(ofl,'   - Beam points Right/Starboard only   ');
3: write(ofl,'   - Beam points Right/Starboard only   ');
4: write(ofl,'   - It is assumed that the beam may be pointed Left or Right  ');
5: write(ofl,'   - Beam can be pointed Left or Right  ');
6: write(ofl,'   - It is assumed that the beam may be pointed Left or Right  ');
7: write(ofl,'   - It is assumed that the beam may be pointed Left or Right  ');
end; {case}
writeln(ofl);
close(xfl);

writeln;
write(ofl,'    Date & Time             Satellite Geocentric             Topocentric Coordinates of Satellite   ');
writeln(ofl,'     Satellite to Target      Pass   ');
write(ofl,'                      Latitude  Longitude  Radius   Node    Azimuth   Elevation  ZenithAngle  Range ');
writeln(ofl,'     Azimuth  NadirAngle   Left/Right');

assign(ifl,'MinimumDistances.txt');
reset(ifl);
assign(pfl,'EphemXYZ.txt');
reset(pfl);
while not(eof(ifl)) do begin
n:=n+1;
if (n mod 20) = 0 then writeln('Tick');
readln(ifl,j,i1,i2,i3,i4,i5,i6,satlat,satlon,satrad,r,c,cc,sdir);
while (k <> j) do readln(pfl,k,yy,mm,day,hh,mn,ss,xs,ys,zs);
dx:=xs-xt;
dy:=ys-yt;
dz:=zs-zt;
dum:=sqrt(sqr(dx) + sqr(dy) + sqr(dz));   {Should equal r}
dlat:=satlat - tarlat;
dlon:=satlon - tarlon;

rs:=(slat*clon*dx) + (slat*slon*dy) - (clat*dz);  {from celetrak.com/columns/v02n02}
re:=-(slon*dx) + (clon*dy);
rz:=(clat*clon*dx) + (clat*slon*dy) + (slat*dz);

elev:=asin(rz/r)*raddeg;
za:=90.0 - elev;
azim:=atan2(-re,rs)*raddeg;
azim:=azim - 180.0;
while (azim < 0.0) do azim:=azim + 360.0;

dum:=sqrt(sqr(rs) + sqr(re) + sqr(rz));   {Should equal r}

sslat:=sin(satlat*degrad);
cslat:=cos(satlat*degrad);
sslon:=sin(satlon*degrad);
cslon:=cos(satlon*degrad);

ps:=(sslat*cslon*dx) + (sslat*sslon*dy) - (cslat*dz);  {from celetrak.com/columns/v02n02}
pe:=-(sslon*dx) + (cslon*dy);
pz:=(cslat*cslon*dx) + (cslat*sslon*dy) + (sslat*dz);

look:=asin(-pz/r)*raddeg;
nadira:=90.0 + look; {look angle away from nadir}
lazim:=atan2(-pe,ps)*raddeg;
while (lazim < 0.0) do lazim:=lazim + 360.0;

dum:=sqrt(sqr(ps) + sqr(pe) + sqr(pz));   {Should equal r}
if (abs(r-dum) > 0.001) then writeln('Error! ',r:12:4,dum:12:4);

shand:=' ';
if ((sdir = 'A') and (lazim < 180.0)) then shand:='R';
if ((sdir = 'A') and (lazim > 180.0)) then shand:='L';
if ((sdir = 'D') and (lazim < 180.0)) then shand:='L';
if ((sdir = 'D') and (lazim > 180.0)) then shand:='R';
if (shand = ' ') then writeln('Error in determining direction of pass/beam');

if ((za > zamin) and (za < zamax)) then begin    {constraints on incident angles = zenith angles}
if ((satno = 1) and (shand = 'R')) then                     {Radarsat-2: may be possible to slew left at some stage}
writeln(ofl,yy:5,mm:3,day:3,hh:3,mn:3,ss:3,satlat:9:2,satlon:10:2,satrad:10:2,'    ',sdir,'  ',
   azim:11:2,elev:10:2,za:11:2,r:11:2,'   ',lazim:10:2,nadira:10:2,'         ',shand,'  ');

if (((satno = 2) or (satno = 3)) and (shand = 'R')) then    {TerraSAR-X and TanDEM-X: points right only}
writeln(ofl,yy:5,mm:3,day:3,hh:3,mn:3,ss:3,satlat:9:2,satlon:10:2,satrad:10:2,'    ',sdir,'  ',
   azim:11:2,elev:10:2,za:11:2,r:11:2,'   ',lazim:10:2,nadira:10:2,'         ',shand,'  ');

if (satno > 3) then
writeln(ofl,yy:5,mm:3,day:3,hh:3,mn:3,ss:3,satlat:9:2,satlon:10:2,satrad:10:2,'    ',sdir,'  ',
   azim:11:2,elev:10:2,za:11:2,r:11:2,'   ',lazim:10:2,nadira:10:2,'         ',shand,'  ');
                                          end; {ZA limits}
                      end; {not eof ifl}

close(ifl);
close(pfl);
close(ofl);

{readln;  }

end.
