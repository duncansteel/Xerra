program tiffstats(input, output, stdErr);

{From Simeon, 2020 July 30}
{Will run in any folder - 2020 October 06}

uses
  fpreadtiff, fpimage;

var
  image: TFPCustomImage;
  reader: TFPCustomImageReader;
  i, j: integer;
  clr: TFPColor;
begin
  if (ParamCount < 1) then
  begin
    writeln('Usage: tiffstats <file.tif>');
    exit;
  end;

  Image := TFPMemoryImage.Create(10, 10);
  Reader := TFPReaderTiff.Create;
  Image.LoadFromFile(paramStr(1), Reader);

    for j := 0 to Image.Height - 1 do
      for i := 0 to Image.Width - 1 do
      begin
        clr := Image.Colors[i, j];
        writeln('i=',i,',j=',j,' red=',clr.red, ' green=',clr.green, ' blue=',clr.blue);
      end;
end.
