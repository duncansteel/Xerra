program EarthRadii(input, output);

{Duncan Steel, 2020 October 06}
{Program to output radius of Earth as a function of latitude for WGS84 ellipsoid}

var
  i : integer;
  pi,degrad,sma,smb,sma2,smb2,r,psi,x,t1,t2,t3,t4 : real;
  ofl : text;

begin
assign(ofl,'EarthRadii.txt');
rewrite(ofl);
writeln; writeln;
pi:=4.0*arctan(1.0);
degrad:=pi/180.0;
sma:=6378137.0;    {Semi-major axis of WGS84 ellipsoid}
smb:=6356752.314245; {Semi-minor axis of WGS84 ellipsoid}
sma2:=sqr(sma);
smb2:=sqr(smb);

for i:=0 to 9000 do begin
x:=i/100.0;
psi:=x*degrad; {latitude, to get radius r}
t1:=sin(psi);
t2:=cos(psi);
t3:=sqr(sma2*t2) + sqr(smb2*t1);
t4:=sqr(sma*t2) + sqr(smb*t1);
r:=sqrt(t3/t4); {Earth radius at latitude psi or i/100}
writeln(ofl,i:6,x:8:2,r:18:6);
                      end; {i=0 to 9000}
writeln(ofl);
writeln('Done - hit return'); writeln;
readln;
close(ofl);

end.
