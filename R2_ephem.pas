program R2_ephem(input, output, stdErr);

{Duncan Steel, 2020 August 07}
{Folders changed 2020 October 06}

{Program for reading in Radarsat-2 ephemeris written out of STK}

var
  i,j,k,m,n : integer;
  r,t,x,y,z,h,v,vx,vy,vz,tstep,t1,t2,dt,tlast,gmearth,rmin,rmax,sma,ecc,v2,v2sqrt : real;
  ifl,ofl,pfl : text;
  sdate : array[1..12] of char;
  c : char;

begin
c:=' ';
for i:=1 to 12 do sdate[i]:=c;
gmearth:=3.9860044188E14; {gravitational constant * mass of the Earth in SI units}

assign(ifl,'Radarsat-2_01-May_1sec.e');
reset(ifl);
assign(ofl,'R2.txt');
rewrite(ofl);
assign(pfl,'R2.csv');
rewrite(pfl);
writeln; writeln;
for i:=1 to 33 do readln(ifl);
for i:=1 to 22 do read(ifl,c);
for i:=1 to 12 do begin read(ifl,c); sdate[i]:=c; end; readln(ifl);
for i:=1 to 12 do write(pfl,sdate[i]); writeln(pfl);
for i:=1 to 12 do write(ofl,sdate[i]);
writeln(ofl); writeln(ofl);
for i:=1 to 3 do readln(ifl);

readln(ifl,t1);
readln(ifl,t2);
tstep:=t2-t1; {Time step in seconds}
reset(ifl);
for i:=1 to 37 do readln(ifl);

t:=0.0;
tlast:=0.0;
rmin:=1.0E12;
rmax:=0.0;

while (t < 86400.0) do begin
readln(ifl,t,x,y,z,vx,vy,vz);
dt:=t - tlast;
if ((abs(dt - tstep) < 0.000001) or (abs(t) < 0.000001)) then begin
tlast:=t;
r:=sqrt(sqr(x) + sqr(y) + sqr(z));
if (r > rmax) then rmax:=r;
if (r < rmin) then rmin:=r;
v:=sqrt(sqr(vx) + sqr(vy) + sqr(vz));
h:=r - (6378137.0);
writeln(ofl,t:10:1,x:12:1,y:12:1,z:12:1,r/1000.0:12:3,h/1000.0:12:3,v/1000.0:12:3);
writeln(pfl,t:10:1,' , ',x:15:3,' , ',y:15:3,' , ',z:15:3,' , ',r/1000.0:15:3,' , ',
            h/1000.0:15:3,' , ',v/1000.0:15:6);
                                     end;
                   end;    {while t < 86400}

sma:=(rmin + rmax)/2.0;
ecc:=1.0 - (rmin/sma);
v2:=gmearth*(1.0/sma);
v2sqrt:=sqrt(v2);
writeln(rmin/1000.0:12:3,rmax/1000.0:12:3,sma/1000.0:12:3,ecc:12:8,v2sqrt:15:6);
writeln; writeln;
readln;
close(ifl);
close(ofl);
close(pfl);

end.
