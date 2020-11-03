program R2tiffstats(input, output, stdErr);

{Duncan Steel, 2020 July 30}
{Folder changed 2020 October 06}
{Program for fiddling with Radarsat-2 data; this is for clipped geotiff image from 2020 May 05, square
covering -37 to -38, 157 to 158}

{Based on routine from Simeon, 2020 July 30}

uses
  fpreadtiff, fpimage;

var
  img : TFPCustomImage;
  rdr : TFPCustomImageReader;
  i,j,k,m,n : integer;
  s,t : real;
  clr : TFPColor;
  x : array[0..500] of real;
  ofl : text;
  answ,c : char;

begin
answ:=' ';
c:=' ';
assign(ofl,'C:\FPC\SARdata\R2tiffstatsoutput.csv');
rewrite(ofl);
for i:=0 to 500 do x[i]:=0.0;
  m:=0;  n:=1000; t:=0.0;
  img := TFPMemoryImage.Create(4598, 4598);
  rdr := TFPReaderTiff.Create;
  img.LoadFromFile('C:\FPC\SARdata\R2_-37-38_157-158.tif',rdr);
  rdr.Free;

     for j := 0 to (img.Height - 1) do
      for i := 0 to (img.Width - 1) do
      begin
        clr := img.Colors[i, j];
        if (clr.red) <> (clr.green) then writeln('Error at ',i:8,j:8);
        if (clr.red) <> (clr.blue) then writeln('Error at ',i:8,j:8);
        if (clr.blue) <> (clr.green) then writeln('Error at ',i:8,j:8);
        k:=clr.red;
        if (k > m) then m:=k;
        if (k < n) then n:=k;
        x[k]:=x[k]+1.0;
        t:=t + k;
      end;

writeln; writeln;
s:=0.0;
for i:=0 to 500 do s:=s+x[i];
{writeln(ofl,'Minimum, maximum, samples, mean = ',n:8,m:8,s:12:1,(t/s):12:1);   }
for i:=0 to 500 do writeln(ofl,i:4,' , ',x[i]:9:1);
writeln('Minimum, maximum, samples, mean = ',n:8,m:8,s:12:1,(t/s):12:1);
writeln; writeln;
readln;
close(ofl);
img.Free;

end.



{
      for j := 0 to 5 do
      for i := 0 to 5 do
      begin
        clr := img.Colors[i, j];
        writeln('i=',i,',j=',j,' red=',clr.red:5, ' green=',clr.green:5, ' blue=',clr.blue:5);
      end;
 }
