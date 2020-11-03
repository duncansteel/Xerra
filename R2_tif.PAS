program R2_tif(input,output);

{Duncan Steel, 2020 July 30}
{Folder changed 2020 October 06}
{Program for fiddling with Radarsat-2 data; this is for clipped geotiff image from 2020 May 05, square
covering -37 to -38, 157 to 158}

var i,j,k,l,m,n : integer;
    a : real;
    x : array[1..2000,1..2000] of integer;
    z : array[0..255] of real;
    ifl : text;
    answ,c : char;

begin
answ:=' ';
a:=0.0;
for i:=0 to 255 do z[i]:=0.0;

assign(ifl,'C:\FPC\SARdata\R2_-37-38_157-158.tif');
reset(ifl);
i:=0;
n:=19166;
for m:=1 to n do read(ifl,c);
l:=n;
{while not(eof(ifl)) do begin }
for m:=1 to 163 do begin
a:=a + 1.0;
i:=i+1;
l:=l+1;
read(ifl,c);
k:=ord(c);
{if k = 26 then writeln(i:6,'  ',k:6);  }
writeln(i:6,l:6,'  ',k:6);
{if (i mod 4) = 0 then begin i:=0; writeln(a:10:1); end; }
                         end;

writeln; writeln;
writeln(a:10:1);
writeln;

i:=0; j:=1;
{
while not eof(ifl) do begin
read(ifl,c);
k:=ord(c);
i:=i+1;

if (i = 2000) then begin
i:=0;
j:=j+1;
                     end;
                      end;
for m:=0 to 20 do begin
write(z[m]:9:1);
if (m mod 10) = 0 then writeln;
                     end;
}
close(ifl);

writeln('Hit return.........  ');
readln;

end. {of program R2_2k}
