procedure ti(var ff:text);
{
 Writes current date and time to file ff. Must declare
 "Uses sysutils;" in main program. File ff must be open
  for writing; "ti" does not write an end-of-line.
}

Var
 ST : TSystemTime; {for SystemTimeToDateTime}
 mmonth:array[1..12] of string[3];

begin {ti}
 mmonth[ 1]:='Jan';
 mmonth[ 2]:='Feb';
 mmonth[ 3]:='Mar';
 mmonth[ 4]:='Apr';
 mmonth[ 5]:='May';
 mmonth[ 6]:='Jun';
 mmonth[ 7]:='Jul';
 mmonth[ 8]:='Aug';
 mmonth[ 9]:='Sep';
 mmonth[10]:='Oct';
 mmonth[11]:='Nov';
 mmonth[12]:='Dec';

 DateTimeToSystemTime(Now,ST);
 With St do
   begin
     write(ff,day:2,' ',mmonth[month],' ',year,hour:3,'h');
     write(ff,minute:3,'m' {,second:3,'s'});
     {may use: millisecond}
   end;
end; {ti}
