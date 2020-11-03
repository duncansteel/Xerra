program R2_fpreadtif(input,output);

{Duncan Steel, 2020 July 30}
{Folder changed 2020 October 06}
{Program for fiddling with Radarsat-2 data; this is for clipped geotiff image from 2020 May 05, square
covering -37 to -38, 157 to 158}

uses fpimage, fpreadtiff;

var i,j,k,l,m,n : integer;
    clr : TFPColor;
    a : real;
    x : array[1..2000,1..2000] of integer;
    ofl : text;
    answ,c : char;
    img : TFPCustomImage;
    rdr : TFPCustomImageReader;

begin
img:=TFPMemoryImage.Create(4598,4598);
rdr:=TFPReaderTiff.Create;
img.LoadFromFile('C:\FPC\SARdata\R2_-37-38_157-158.tif',rdr);
rdr.Free;

writeln; writeln;
answ:=' ';
c:=' ';
a:=0.0;

assign(ofl,'C:\FPC\SARdata\R2_fpreadtiffoutput.txt');
rewrite(ofl);

clr:=img.Colors[1,1];

m:=img.Height;
n:=img.Width;

writeln(m:10,n:10);
writeln(ofl,m:10,n:10);

k:=10;
{for i:=1 to k do write(ofl,clr);}

close(ofl);
img.Free;

writeln; writeln;
writeln('Hit return.........  ');
readln;

end. {of program R2_fpreadtif}
