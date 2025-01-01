 procedure readstr(var ff:text; var ssout:string);
{
 To read a string from file ff. String ends when
 blank or tab or end-of-line is encountered.
 Skips white space (blank or tab) preceding string.
 Warning: sysutils also contains a procedure,
 readstr, with a different function.
}
 label 33,66;
 var
  ii:integer;
  ccc:char;
  first:boolean;
 begin
  first:=true; ii:=0;

  repeat
   inc(ii);
   33:
   if eoln(ff) {do not use seekeoln!}
   then if first
    then begin ii:=1; goto 66; end
    else goto 66
   else {not eoln} read(ff,ccc);
   {Skip white space preceding string}
   if first then if (ccc=' ') or (ccc=tab) then goto 33;
   first:=false;
   ssout[ii]:=ccc;
  until (ssout[ii]=' ') or (ssout[ii]=tab);

  66:
  {ssout[0]:=chr(ii-1);} {length = ii-1, not counting trailing blank}
  setlength(ssout,ii-1);
 end; {readstr.p}
