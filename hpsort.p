
procedure hpsort(n:integer; var ra:vector);
{
 Sorts an array ra(1:n) into ascending numerical
 order using the Heapsort algorithm. n is input;
 ra (real) is replaced on output by its sorted
 rearrangement. Based on Press et al, "Numerical
 Recipes in F77", 1992, p. 329

 On return, the smallest number will be ra[1],
 the largest ra[n]; values in increasing order.
}

label 10,20,95;
var
 i,ir,j,L:integer;
 rra:real;

begin {hpsort}
 if n<2 then goto 95;
 L := 1 + (n div 2);
 ir := n;

 10:
 if (L>1) then begin
  L:=L-1;
  rra:=ra[L];
 end else begin
  rra:=ra[ir];
  ra[ir]:=ra[1];
  ir:=ir-1;
  if (ir=1) then begin
   ra[1]:=rra;
   goto 95;
  end;
 end;
 i:=L;
 j:=L+L;

 20:
 if (j <= ir) then begin
  if (j<ir) then begin
   if (ra[j]<ra[j+1]) then inc(j);
  end;
  if (rra<ra[j]) then begin
   ra[i]:=ra[j];
   i:=j;
   j:=j+j;
  end else j:=ir+1;
  goto 20;
 end; {if}

 ra[i]:=rra;
 goto 10;
 95:
end; {hpsort}
