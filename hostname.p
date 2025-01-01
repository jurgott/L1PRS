procedure PCname(var ff:text);
{
 Writes the hostname to file ff. Result will be, for example:
 COMPUTERNAME=64GB-PC
 Empty code for Windows
}
const
 sproc='proc';

Var
 ftmp:text;
 sss:string;
 s4:string[4];
 nCPU:integer;
 
begin
{$ifdef unix}
 assign(ftmp,'/etc/hostname'); reset(ftmp);
 read(ftmp,sss);
 close(ftmp);
 writeln(ff,'Computername = ',sss);
{Number of CPUs, Linux only}
 nCPU:=0;
 assign(ftmp,'/proc/cpuinfo'); reset(ftmp);
 while not eof(ftmp) do begin
  readln(ftmp,s4);
  if s4=sproc then inc(nCPU);
 end;
 close(ftmp);
 writeln(ff,'   CPU count = ',nCPU);
{$endif}
end; {PCname}
