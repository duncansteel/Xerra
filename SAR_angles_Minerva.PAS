program SAR_angles_Minerva(input,output);

{Duncan Steel, 2020 October 13}
{Program for deriving incident angles for SAR from different satellites given a range of look angles}
{This is for the latitude of the Minerva reefs: -23 deg}

var i,k : integer;
    pi,convrd,h,r,x,a1,a2,b1,b2,c1,c2,looka,ranga : real;
    ofl : text;

{------------------------------------------------------------------}

FUNCTION ASIN(XYZ : REAL) : REAL;

{Returns the inverse sine of parameter XYZ}

BEGIN
IF ABS(XYZ) > 0.999999 THEN ASIN:=XYZ*PI/2.0
     ELSE ASIN:=ARCTAN(XYZ/SQRT(1.0-(XYZ*XYZ)));
END;   {of function ASIN}

{------------------------------------------------------------------}

begin {main program}
pi:=4.0 * arctan(1.0); convrd:=180.0/pi;
Assign(ofl,'SAR_angles.txt');
rewrite(ofl);

r:=6375.0; {Assumed Earth radius, appropriate for latitude 22 degrees}
writeln; writeln;
write('Assumed altitude (km)?  '); readln(h);
{At latitude 40 degrees: }
{h = 712 km for Sentinel-1A}  {Varies by about 1 km due to finite e}
{h = 701 km for Sentinel-1B}  {Varies by about 1 km due to finite e}
{h = 525 km for TerrSAR-X}    {Varies by about 1 km due to finite e}
{h = 808 km for Radarsat-2}   {Varies by about 1 km due to finite e}
{h = 603 km for NovaSAR-1}    {Varies by about 6 km due to finite e}

writeln;
x:=r / (r + h);
write('Smaller incident angle? '); readln(a1);
writeln;
write('Larger incident angle?  '); readln(a2);
writeln;
c1:=180.0 - a1;
c2:=180.0 - a2;
b1:=(asin(x * sin(c1/convrd)))*convrd;
b2:=(asin(x * sin(c2/convrd)))*convrd;
writeln;
writeln('Beam limits',(90.0-b2):6:2,(90.0-b1):6:2);
looka:=(b1 + b2)/2.0;
ranga:=(b2 - b1)/2.0;

writeln;
writeln('Look angle =  ',(90.0-looka):6:2,' +/- ',ranga:6:2);
writeln(ofl,'Look angle =  ',(90.0-looka):6:1,' +/- ',ranga:6:1);
writeln;
readln;

close(ofl);

end.  {main program}
