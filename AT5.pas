program AT5(input, output, stdErr);

{Duncan Steel, 2020 October 30}

{AT = Auto Tasking}

{Based on R2_LLR.pas...}
{Program for reading in Radarsat-2 ephemeris written out of STK, this time in Fixed LLR coordinates}
{Note latitude is GEOCENTRIC not GEODETIC}   {In STK, LLA implies geodetic}

var
  j,k : LongInt;
  i,n,day,mm,yy,hh,mn,ss,i1,i2,i3,i4,i5,i6 : integer;
  r,ctarlat,tarlat,tarlon,xt,yt,zt,xs,ys,zs,slat,clat,sma,smb,gmearth,e2,fss,lon,rlast,
  satrad,pi,degrad,raddeg,sma2,smb2,tarrad,satlat,satlon,b2oa2,dum,t1,t2,t3,t4,dx,dy,dz : real;
  ifl,ofl,pfl : text;
  m : array[1..3] of char;
  s : array[1..50] of char;
  c : char;

begin
c:=' ';
for i:=1 to 3 do m[i]:=c;
for i:=1 to 50 do s[i]:=c;
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
tarrad:=tarrad/1000.0; {in km}

assign(ifl,'EphemLLR.txt');
reset(ifl);
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
write(ofl,'Target geocentric lat, lon, radius:',ctarlat*raddeg:9:4,tarlon*raddeg:10:4,tarrad:10:4);
writeln(ofl,'; X,Y,Z: ',xt/1000.0:14:6,yt/1000.0:14:6,zt/1000.0:14:6);
for n:=1 to (i-2) do write(ofl,s[n]);
writeln(ofl);
write(ofl,'  Index    Year   Mon   Day  Hour   Min   Sec   ');
writeln(ofl,'Geocentric Lat         Lon         Distance     Target Distance');

j:=-1; {counter}
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
while not(eof(ifl)) do begin
readln(ifl,j,yy,mm,day,hh,mn,ss,satlat,satlon,satrad,r);
if (j > 0) then begin
if (r < 3000.0) then begin
if (r < rlast) then begin
rlast:=r;
                      end {r < rlast}
                      else begin
                      if (n = 0) then
                      begin writeln(ofl,(j-1):6,i1:6,i2:4,i3:4,i4:4,i5:4,i6:4,t1:14:6,t2:14:6,t3:14:6,t4:14:6);
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
for i:=1 to 74 do read(ifl,c);
read(ifl,xt,yt,zt);
close(ifl);

j:=0;
k:=-1;
n:=0;
assign(ifl,'MinimumDistances.txt');
reset(ifl);
assign(pfl,'EphemXYZ.txt');
reset(pfl);
assign(ofl,'OutPut.txt');
rewrite(ofl);

while not(eof(ifl)) do begin
n:=n+1;
readln(ifl,j,i1,i2,i3,i4,i5,i6,satlat,satlon,satrad,r);
while (k <> j) do readln(pfl,k,yy,mm,day,hh,mn,ss,xs,ys,zs);
dx:=xs-xt;
dy:=ys-yt;
dz:=zs-zt;
dum:=sqrt(sqr(dx) + sqr(dy) + sqr(dz));   {Should equal r}
writeln(ofl,n:3,yy:6,mm:4,day:4,hh:4,mn:4,ss:4,
     xt:14:6,yt:14:6,zt:14:6,xs:14:6,ys:14:6,zs:14:6,r:14:6,satlat:14:6,satlon:14:6,satrad:14:6);
                         end;

close(ifl);
close(pfl);
close(ofl);

{readln;   }

end.



writeln(n,r:14:6,dum:14:6);