program Scale(input,output);

{Duncan Steel, 2014 January 02}
{Runs in any folder: 2020 October 06}

{Makes a BMP scale}

{Starting code may be useful in making other bitmaps}

var i,j,k,m,n : integer;
    o1,o2,o3 : array[1..256,1..10] of char;
    ofl : text;

begin

assign(ofl,'Scale.bmp');
rewrite(ofl);
write(ofl,char(66),char(77));  {BM for Windows}
write(ofl,char(54),char(10),char(0),char(0)); {File size 2614 bytes: first
   is out of 256 bytes, then * 256, then *256*256, final is *256*256*256 }
write(ofl,char(0),char(0),char(0),char(0));
write(ofl,char(54),char(0),char(0),char(0));
   {Total header is 54 bytes including above 14 bytes}
write(ofl,char(40),char(0),char(0),char(0)); {Header is 40 bytes}
write(ofl,char(0),char(1),char(0),char(0)); {Image width 256 pixels}
write(ofl,char(10),char(0),char(0),char(0)); {Image height PLUS 10 pixels}
write(ofl,char(1),char(0),char(24),char(0)); {One colour pane; 24 bits per pixel}
write(ofl,char(0),char(0),char(0),char(0)); {No compression}
write(ofl,char(0),char(10),char(0),char(0)); {Image size: 2560 pixels here}
write(ofl,char(100),char(0),char(0),char(0)); {Horizontal resolution: 100 pixels per metre here}
write(ofl,char(100),char(0),char(0),char(0)); {Vertical resolution: 100 pixels per metre here}
write(ofl,char(0),char(0),char(0),char(0)); {0 defaults to 2^n colours in palette}
write(ofl,char(0),char(0),char(0),char(0)); {0 defaults to all colours are important}
{Above makes BMP 256 x 10, each pixel 10mm}

for k:=1 to 10 do begin
for i:=1 to 256 do begin
m:=i-1;

if (m <= 127) then n:=255 else n:=256-(2*(m-127));
if (n < 0) or (n > 255) then writeln(i:6,k:6,n:6);
o1[i,k]:=char(n);

if (m <= 127) then n:=(2*m)+1 else n:=(2*(255-m))+1;
if (n < 0) or (n > 255) then writeln(i:6,k:6,n:6);
o2[i,k]:=char(n);

if (m <= 127) then n:=2*m else n:=255;
if (n < 0) or (n > 255) then writeln(i:6,k:6,n:6);
o3[i,k]:=char(n);

end; {j=1 to 256}
end; {k=1 to 10}

for j:=1 to 2 do begin  {Tick marks}
for k:=1 to 17 do begin
m:=((k-1) * 15) + 8;         {This is for a maximum of 8 W/m2}
o1[m,j]:=char(0);
o2[m,j]:=char(0);
o3[m,j]:=char(0);
                   end;
                   end;

for k:=1 to 10 do begin
for i:=1 to 256 do begin
write(ofl,o1[i,k],o2[i,k],o3[i,k]);
end; {j=1 to 256}
end; {k=1 to 10}

close(ofl);

writeln('Hit return.........  ');
readln;

end. {of program Scale.pas}
