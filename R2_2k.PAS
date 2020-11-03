program R2_2k(input,output);

{Duncan Steel, 2020 July 29}
{Folder changed 2020 October 06}
{Program for fiddling with Radarsat-2 data; this is for 2k x 2k clipped image from 2020 May 05}

var i,j,k,m,n : integer;
    x : array[1..2000,1..2000] of integer;
    z : array[0..255] of real;
    ifl : text;
    answ,c : char;

begin
answ:=' ';

for i:=0 to 255 do z[i]:=0.0;

assign(ifl,'C:\FPC\SARdata\R2_2k.bmp');
reset(ifl);
for i:=1 to 118 do begin
read(ifl,c); {Header}
{write(ord(c):6);
if ((i+2) mod 4) = 0 then writeln;}
end;

i:=0; j:=1;
while not eof(ifl) do begin
read(ifl,c);
k:=ord(c);
i:=i+1;
if (k = 0) then begin x[i,j]:=0; x[i+1,j]:=0; end;
if (k = 1) then begin x[i,j]:=0; x[i+1,j]:=1; end;
if (k = 16) then begin x[i,j]:=1; x[i+1,j]:=0; end;
if (k = 17) then begin x[i,j]:=1; x[i+1,j]:=1; end;
i:=i+1;
if (k = 0) then z[0]:=z[0] + 2.0;
if (k = 1) then begin z[0]:=z[0] + 1.0; z[1]:=z[1] + 1.0; end;
if (k = 16) then begin z[0]:=z[0] + 1.0; z[1]:=z[1] + 1.0; end;
if (k = 17) then z[1]:=z[1] + 2.0;

if (i = 2000) then begin
{if (j mod 40) = 0 then readln;}
{writeln(i:5,j:5);}
i:=0;
j:=j+1;
end;
                    end;
for m:=0 to 20 do begin
write(z[m]:9:1);
if (m mod 10) = 0 then writeln;
                     end;
close(ifl);
writeln('Hit return.........  ');
readln;

end. {of program R2_2k}
