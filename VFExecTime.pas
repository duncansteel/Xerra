program VFExecTime(input, output);

{Duncan Steel, 2020 October 06}
{Output folder changed 2020 October 06}

{Program for solving the Vincenty forward/direct problem}
{This version for testing the execution time}

uses DateUtils,SysUtils;
var
  t00,t01 : TDateTime;
  i,j,k,n : integer;
  pi,piby2,twopi,degrad,raddeg,gmearth,fo16,sma,smb,sma2,smb2,dsm2,e2,oneof,onemf,aob,f,dist,alpha1,alpha2,a2b2b2,
  la1,la2,lo1,lo2,dstep,lab,lob,alphab,avt : real;
  ofl,dfl : text;
  c : char;

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
procedure vincentyinverse(lat1,lon1,lat2,lon2 : real; var alph1,alph2,s : real);
var u1,u2,bigl,lam,alpha,sigma,sig1,sigm,q1,q2,s1,s2,c1,c2,dlam,lamlast,ssig,csig,su1,cu1,su2,cu2,deltsig,
    cu1su2,su1cu2,su1su2,cu1cu2,c2sigm,slam,clam,salph,csqalph,bigc,biga,bigb,usq,sqc2sigm,k1 : real;
begin
u1:=1000.0;     {Check values}
u2:=1000.0;
q1:=lat1*degrad;
q2:=lat2*degrad;
s1:=sin(q1);
s2:=sin(q2);
c1:=cos(q1);  {latitude is between -90 and +90 so cosine is always positive}
c2:=cos(q2);
if (c1 > 0.0) then u1:=arctan(onemf*s1/c1) else begin
                                                if (s1<-0.999999) then u1:=-piby2;  {South pole}
                                                if (s1>0.999999) then u1:=-piby2;   {North pole}
                                                  end;
if (u1 > 999.9) then writeln('u1 error in Vincenty Inverse');
su1:=sin(u1);
cu1:=cos(u1);
if (c2 > 0.0) then u2:=arctan(onemf*s2/c2) else begin
                                                if (s2<-0.999999) then u2:=-piby2;  {South pole}
                                                if (s2>0.999999) then u2:=-piby2;   {North pole}
                                                  end;
if (u2 > 999.9) then writeln('u2 error in Vincenty Inverse');
su2:=sin(u2);
cu2:=cos(u2);
cu1su2:=cu1*su2;
su1cu2:=su1*cu2;
su1su2:=su1*su2;
cu1cu2:=cu1*cu2;

bigl:=(lon1 - lon2)*degrad;
dlam:=1000.0;
lamlast:=1000.0;
lam:=bigl;
while (abs(lam-lamlast) > 1.0e-12) do begin
lamlast:=lam;
slam:=sin(lam);
clam:=cos(lam);
ssig:=sqrt((sqr(cu2*slam)) + (sqr(cu1su2 - (su1cu2*clam))));
csig:=su1su2 + (cu1cu2*clam);
sigma:=atan2(ssig,csig);
salph:=cu1cu2*slam/ssig;
csqalph:=1.0 - sqr(salph);
c2sigm:=csig - (2.0*su1su2/csqalph);
sqc2sigm:=sqr(c2sigm);
bigc:=(fo16*csqalph)*(4.0 + (f*(4.0 - (3.0*csqalph))));
lam:=bigl + ((1.0-bigc)*f*salph*(sigma + (bigc*ssig*(c2sigm + (bigc*csig*(-1.0 + (2.0*sqr(c2sigm))))))));
                                       end; {while... converges}
slam:=sin(lam);
clam:=cos(lam);
usq:=csqalph*a2b2b2;

{Original Vincenty formulae for A & B:
biga:=1.0 + ((usq/16384.0) * (4096.0 + (usq * ((usq * (320.0 - (175.0 * usq))) - 768.0))));
bigb:=(usq/1024.0)*(256.0 + (usq*((usq*(74.0-(47.0*usq))) - 128.0)));
}

{Vincenty's modification to derive A & B: }
q1:=sqrt(1.0 + usq); {re-use variable}
k1:=(q1 - 1.0)/(q1 + 1.0);
q2:=0.25*sqr(k1); {re-use variable}
biga:=(1.0 + q2) / (1.0 - k1);
bigb:=k1 * (1.0 - (1.5*q2));

deltsig:=bigb*ssig*(c2sigm + ((bigb/4.0)*((csig*(((2.0*sqc2sigm)-1.0)) -
                  ((bigb/6.0)*c2sigm*((4.0*sqr(ssig))-3.0)*((4.0*sqc2sigm)-3.0))))));
s:=smb*biga*(sigma - deltsig);
q1:=cu2*slam;               {re-use variables}
q2:=cu1su2 - (su1cu2*clam);
alph1:=atan2(q1,q2);
alph1:=twopi - alph1;
while (alph1 > twopi) do alph1:=alph1-twopi;
while (alph1 < 0.0) do alph1:=alph1+twopi;
q1:=cu1*slam;               {re-use variables}
q2:=-su1cu2 + (cu1su2*clam);
alph2:=atan2(q1,q2);
alph2:=twopi - alph2;
while (alph2 > twopi) do alph2:=alph2-twopi;
while (alph2 < 0.0) do alph2:=alph2+twopi;
alph1:=alph1*raddeg; {in degrees}
alph2:=alph2*raddeg;

end; {of procedure vincentyinverse}

{------------------------------------------------------------------}

procedure vincentyforward(lat1,lon1,alph1,pstep : real; var lan,lon,alphn : real);

var u1,tu1,salph1,bigl,lam,alpha,sigma,sig1,twosigm,ctwosigm,c2twosigm,q1,q2,sqsalph,calph1,d1,d2,d3,d4,
    s1,c1,ssig,csig,su1,cu1,deltsig,salph,csqalph,bigc,biga,bigb,usq,k1,psa,b4,b6,dum,slast : real;
    nn : integer;

begin
u1:=1000.0;     {Check values}
alph1:=alph1*degrad; {in radians}
q1:=lat1*degrad;
s1:=sin(q1);
c1:=cos(q1);  {latitude is between -90 and +90 so cosine is always positive}
if (c1 > 0.0) then begin
                   tu1:=onemf*s1/c1;
                   u1:=arctan(tu1);
                   su1:=sin(u1);
                   cu1:=cos(u1);
                     end
                     else begin
                     if (s1<-0.999999) then begin u1:=-piby2; su1:=-1.0; cu1:=0.0; tu1:=-1.0e10; end; {South pole}
                     if (s1>0.999999) then begin u1:=piby2; su1:=1.0; cu1:=0.0; tu1:=1.0e10; end;   {North pole}
                            end;
if (u1 > 999.9) then writeln('u1 error in Vincenty Forward');
salph1:=sin(alph1);
calph1:=cos(alph1);
sig1:=atan2(tu1,calph1);
salph:=cu1 * salph1;
sqsalph:=sqr(salph);
csqalph:=1.0 - sqsalph;
usq:=(1.0 - sqsalph)*a2b2b2;
{Vincenty's modification to derive A & B: }
q1:=sqrt(1.0 + usq); {re-use variable}
k1:=(q1 - 1.0)/(q1 + 1.0);
q2:=0.25*sqr(k1);
biga:=(1.0 + q2) / (1.0 - k1);
bigb:=k1 * (1.0 - (1.5*q2));
b4:=bigb/4.0;
b6:=bigb/6.0;
psa:=pstep / (smb*biga); {step s divided by smb times A}
sigma:=psa;      {initial guess for sigma}
deltsig:=1000.0; {large start value}
slast:=1000.0; {large start value}
nn:=0; {Counter}
while ((abs(sigma-slast) > 1.0E-12) and (nn < 11)) do begin
nn:=nn+1;
twosigm:=(2.0*sig1) + sigma;
ctwosigm:=cos(twosigm);
c2twosigm:=sqr(ctwosigm);
dum:=2.0*c2twosigm;
ssig:=sin(sigma);
csig:=cos(sigma);
d1:=bigb*ssig;
d2:=csig*(dum - 1.0);
d3:=b6 * ctwosigm * ((4.0*sqr(ssig)) - 3.0) * ((2.0*dum) - 3.0);
d4:=b4 * (d2 - d3);
deltsig:=d1 * (ctwosigm + d4);
slast:=sigma;
sigma:=psa + deltsig;
                                end; {while}
ssig:=sin(sigma);
csig:=cos(sigma);
d1:=(su1*csig) + (cu1*ssig*calph1);
d2:=sqr((su1*ssig) - (cu1*csig*calph1));
d3:=sqrt(sqsalph + d2);
d4:=onemf*d3;
d2:=atan2(d1,d4);
lan:=d2*raddeg;  {in degrees}
d1:=ssig*salph1;
d2:=(cu1*csig) - (su1*ssig*calph1);
lam:=atan2(d1,d2);
bigc:=fo16*csqalph*(4.0 + (f*(4.0 - (3.0*csqalph))));
d4:=bigc*csig*(dum - 1.0);
d3:=ctwosigm + d4;
d2:=bigc * ssig * d3;
d1:=sigma + d2;
d4:=(1.0 - bigc) * f * salph * d1;
bigl:=lam - d4;
lon:=(bigl*raddeg) + lon1;  {in degrees}
while (lon < 0.0) do lon:=lon + 360.0;
while (lon > 360.0) do lon:=lon - 360.0;

d1:=(cu1*csig*calph1) - (su1*ssig);
d2:=atan2(salph,d1);
alphn:=d2;
while (alphn > twopi) do alphn:=alphn-twopi;
while (alphn < 0.0) do alphn:=alphn+twopi;
alphn:=alphn*raddeg; {in degrees}

end; {of procedure vincentyforward}

{------------------------------------------------------------------}

begin
assign(dfl,'VFExecTime.txt');
rewrite(dfl);

avt:=0.0;
for k:=1 to 1000 do begin
t00:=Now;
assign(ofl,'Dummy.txt');
rewrite(ofl);
writeln; writeln;

c:=' ';
pi:=4.0*arctan(1.0);
piby2:=pi/2.0;
twopi:=2.0*pi;
degrad:=pi/180.0;
raddeg:=180.0/pi;
gmearth:=3.9860044188E14; {gravitational constant * mass of the Earth in SI units}
sma:=6378137.0;    {Semi-major axis of WGS84 ellipsoid}
smb:=6356752.314245; {Semi-minor axis of WGS84 ellipsoid}
aob:=sma/smb;
oneof:=298.257223563; {Reciprocal of WGS84 ellipsoid flattening}
f:=1.0/oneof;
fo16:=f/16.0;
e2:=1.0 - sqr(smb/sma); {Square of eccentricity of WGS84 ellipsoid}
sma2:=sqr(sma);
smb2:=sqr(smb);
dsm2:=sma2 - smb2;
a2b2b2:=dsm2/smb2;
onemf:=1.0-f;

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

vincentyinverse(la1,lo1,la2,lo2,alpha1,alpha2,dist);        {This renders the initial azimuth, alpha1}
writeln(ofl,'(Vincenty Inverse Problem solution)');
writeln(ofl,'Overall: Start latitude Start longitude Start azimuth  End latitude   End longitude   End azimuth    Distance');
writeln(ofl,'          ',la1:7:2,'        ',lo1:7:2,'       ',alpha1:12:6,
     '    ',la2:7:2,'        ',lo2:7:2,'       ',alpha2:12:6,(dist/1000.0):12:3);
writeln(ofl);
i:=1;

n:=100;
dstep:=dist/n; {Each step is 1/n of the total distance, in metres}
writeln(ofl,'(Vincenty Direct Problem solution)');
writeln(ofl,'With small steps of length ',dstep/1000.0:10:6,' km:');
writeln(ofl,'  Step   Start latitude Start longitude Start azimuth  End latitude   End longitude   End azimuth');

for i:=1 to n do begin
vincentyforward(la1,lo1,alpha1,dstep,lab,lob,alphab);       {for a small step along the path}
writeln(ofl,i:6,la1:15:6,lo1:15:6,alpha1:15:6,lab:15:6,lob:15:6,alphab:15:6);
la1:=lab;
lo1:=lob;
alpha1:=alphab;
                    end; {i=1 to n}
writeln; writeln;
writeln(ofl);
close(ofl);

t01:=Now;
writeln(dfl,k:6,'    ',(t01-t00):12);
avt:=avt + (t01-t00);
   end; {k=1 to 1000}
writeln(dfl,'Average execution time = ',(avt/1000.0):12);
close(dfl);

end.
