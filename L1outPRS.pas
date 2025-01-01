{ $ R + }
program L1outPRS(input,output);
{
 (C) 2021-2025 Jurg Ott   GNU General Public License v3

 LINUX: Do cross-validation by the leave-one-out method
 (Agresti 2019) for Polygenic Risk Score (PRS) of single
 variants.
   Principle: Individual is called a "case" if his score
 is above the 95th percentile of control scores.
   If plink issues notes saying that a line in a score file
 was omitted because of allele incompatiblity, try using the
 same dataset but with increased MAF. Even though the full
 dataset may not contain monomorphic SNPs, removing an individual
 may make frequency of allele 1 zero, but this is likely to
 occur only in very small datasets.

 Sort variable is (1-p) [1] or OR [2].
 Score variable is ACC [1] or PPV [2].

20 Mar 2024 written
25 Mar 2024 Use sorted pattern file
26 Mar 2024 shpat = input
10 Apr 2024 new layout
12 Apr 2024 Allow Linux
13 Apr 2024 Work with fields, not columns
14 Apr 2024 append to outfile
15 Apr 2024 Minor allele = reference allele, A1=1
18 Apr 2024 Output changes
28 Apr 2024 Set bounds on pp; see "28 Apr 2024"
30 Apr 2024 "L1outPRS" prefix for all output files
15 Jun 2024 Directly input metric variable
18 Jun 2024 add number of SNPs in output
20 Jun 2024 more detailed output to logfile
22 Jun 2024 add comments on command line
15 Jul 2024 major overhaul
22 Jul 2024 do only one, ll=1, 95th percentile
25 Jul 2024 add hostname, random values of PPV and ACC
26 Jul 2024 minor changes
05 Aug 2024 compute random ACC correctly
10 Dec 2024 syntax error corrected
11 Dec 2024 minor changes
12 Dec 2024 score variable = ACC instead of PPV
15 Dec 2024 make unified program
27 Dec 2024 Rename scoreOR => scorePPV
28 Dec 2024 Fix path to plink19
29 Dec 2024 Test for existence of plink19
01 Jan 2025 Adjust program constants
}
uses sysutils;

type
 integer=longint; real=double;

const
 version='01 Jan 2025';
 maxind=10000;
 maxsnp=13; {change this for different numbers of Nsnp groups}
 two=2.0; bl=' '; one=1.0; half=0.5;
 ndec=3;   {number of decimals in output}
 ll=1; {number of limits for percentiles}
 progname='L1outPRS'; plinkname='/usr/local/bin/plink19';

type
 vector=array[1..maxind] of real;

var
 outfile,reportfile,logfile:text;
 outfilename,reportfilename,logfilename,fnameplink,namescore,
  namesort,fnameplinkbin:string;
 tab{,bs}:char;
 exitcode,iSNP,NNsnp,Nindall,leave,casealleles,ctrlalleles:integer;
 Nsnp:array[1..maxsnp] of integer;
 phenoOrig:array[1..maxind] of shortint;
 scoreOrig:vector;
 FID,IID:array[1..maxind] of string;
 percent:array[1..ll] of real; {[0.90]  0.95  0.99}
 {                Nsnp  %ile phen pred}
 table:array[1..maxsnp,1..ll,1..2,1..2] of integer;
 sortbyOR,scorePPV:boolean;
 pheno:array[1..2] of integer;


{$I readstr.p}
{$I ti.p}


 procedure init;
 var
  i1,i2,i3,i4:integer;
  ok:boolean;

{$I hostname.p}	{calls PCname procedure}

 begin {init}
  writeln;
  writeln('Program ',progname,' version ',version);
  writeln(progname,' performs cross-validation by the');
  writeln('leave-one-out method for GWAS data as implemented');
  writeln('in plink with the --score function.');
  writeln;
  writeln('Plink version 1.9 must reside as an executable');
  writeln('"plink19" in the /usr/local/bin folder');
  writeln;
  writeln('Maximum program constants:');
  writeln(maxind:8,' individuals');
  writeln;
  writeln('Run this program with 3 command line parameters:');
  writeln(' 1) The plink fileset name (without ".map", ".ped")');
  writeln(' 2) Sorting of variants by p-value [1] or odds ratio [2]');
  writeln(' 3) score variable = ACC [1] or PPV [2]');

  if not FileExists(plinkname) then begin
   writeln('Program "plink19" not found in /usr/local/bin');
   halt;
  end;

  if paramcount=0 then halt;
  if (paramcount < 3) then begin
   writeln('There must be 3 command line parameters');
   halt;
  end;
  tab:=chr(9); {bs:=chr(8);}
  fnameplink:=paramstr(1);

  ok := (paramstr(2)='1') or (paramstr(2)='2');
  if not ok then begin
   writeln('Second parameter must be 1 or 2');
   halt;
  end;

  ok := (paramstr(3)='1') or (paramstr(3)='2');
  if not ok then begin
   writeln('Third parameter must be 1 or 2');
   halt;
  end;

  namesort  := paramstr(2);
  namescore := paramstr(3);
  sortbyOR  := (namesort='2');  {sort by OR if true, by (1-p) if false}
  scorePPV  := (namescore='2'); {score = PPV if true, = ACC if false}

  ok := FileExists(fnameplink+'.ped') and FileExists(fnameplink+'.map');
  if not ok then begin
   writeln(fnameplink+'.ped',' or ',fnameplink+'.map',' not found');
   halt;
  end;

  reportfilename := progname+'-'+fnameplink+'-'+namesort+namescore+'.rpt';
  assign(reportfile,'../'+reportfilename); rewrite(reportfile);
  writeln(reportfile,progname+' version ',version,': Polygenic Risk Score for "',fnameplink,'-',namesort,namescore,'" data');
  close(reportfile);

  {            Nsnp      %ile phen pred}
  {table:array[1..maxsnp,1..3,1..2,1..2] of integer;}
  for i1:=1 to maxsnp do {maxsnp=13}
  for i2:=1 to ll do
  for i3:=1 to 2 do
  for i4:=1 to 2 do table[i1,i2,i3,i4]:=0;

  Nsnp[ 1]:=    5;
  Nsnp[ 2]:=   10;
  Nsnp[ 3]:=   20;
  Nsnp[ 4]:=   50;
  Nsnp[ 5]:=  100;
  Nsnp[ 6]:=  200;
  Nsnp[ 7]:=  500;
  Nsnp[ 8]:= 1000;
  Nsnp[ 9]:= 2000;
  Nsnp[10]:= 5000;
  Nsnp[11]:=10000;
  Nsnp[12]:=20000;
  Nsnp[13]:=50000;

  percent[1]:=0.95;
{
  percent[1]:=0.95;
  percent[2]:=0.99;
}
  outfilename:=progname+'.'+fnameplink+'-'+namesort+namescore+'.out';
  logfilename:=progname+'.'+fnameplink+'-'+namesort+namescore+'.log';
  assign(outfile,outfilename);
  rewrite(outfile);
  assign(logfile,logfilename);
  rewrite(logfile);
  writeln(logfile,'Program ',progname,' version ',version);
  writeln(logfile,progname,' does cross-validation by the leave-one-out method');
  writeln(logfile,'for GWAS data as implemented in plink''s --score function.');
  writeln(logfile);
  writeln(logfile,'Maximum program constants:');
  writeln(logfile,maxind:11,' individuals');
  writeln(logfile);
  writeln(logfile,'Input data: ',fnameplink);

  if sortbyOR
  then writeln(logfile,' Sort variable = odds ratio')
  else writeln(logfile,' Sort variable = (1-p)-value');

  if scorePPV
  then writeln(logfile,'Score variable = positive predictive value, PPV')
  else writeln(logfile,'Score variable = prediction accuracy, ACC');

  PCname(logfile);
  ti(logfile); writeln(logfile,': Start');

  for i1:=1 to ll do
  for i2:=1 to 7 do write(outfile,tab,percent[i1]:4:2);
  writeln(outfile);

  write(outfile,'Nsnp');
  for i1:=1 to ll do
  write(outfile,tab,'a',tab,'b',tab,'c',tab,'d',tab,'PPV',tab,'NPV',tab,'ACC');
  writeln(outfile);
  close(outfile);
 end; {init}


{$I hpsort.p} {needed for distributions of scores in cases and controls}



 procedure makeLNor;
{
 Check values of score
 Write assocLNorig.txt:	SNP  A1  PPV
}
 var
  f1,f2:text;
  ss4,sss,A1:string;
  scovar:real; {ACC}
  allsnps,ii,aa,bb,cc,dd,N1:integer; {N1 = number of A1 alleles}
  ok:boolean;
  pp,odds:real;

 begin {makeLNor}

  {Make binary plink files}
  fnameplinkbin:='ppp';
  ss4:=' --file '+fnameplink+' --make-bed --out '+fnameplinkbin;
  exitcode:=ExecuteProcess(plinkname,ss4);
  if exitcode<>0 then begin
   writeln('Exit code = ',exitcode,': Run plink --make-bed');
   halt;
  end;

  {Count SNPs in mapfile}
  assign(f1,fnameplink+'.map'); reset(f1);
  allsnps:=0;
  while not eof(f1) do begin
   inc(allsnps);
   readln(f1);
  end;
  close(f1);
  writeln(logfile,allsnps,' variants in file ',fnameplink+'.map');

  {Count pheno in pedfile}
  for ii:=1 to 2 do pheno[ii]:=0;
  assign(f1,fnameplink+'.ped'); reset(f1);

  while not eof(f1) do begin
   for ii:=1 to 5 do readstr(f1,sss); {FID IID PID MID sex; skip}
   readln(f1,ii); {pheno}
   ok := (ii=1) or (ii=2);
   if not ok then begin
    writeln('Phenotype other than 1 or 2 encountered');
    halt;
   end;
   inc(pheno[ii]);
  end;

  close(f1);
  pp:=pheno[2]/(pheno[1]+pheno[2]);
  writeln(logfile,pheno[2],' cases, ',pheno[1],' controls');
  writeln(logfile,'Proportion of cases = ',pp:6:4);
  casealleles:=2*pheno[2]; ctrlalleles:=2*pheno[1];

  {Run initial assoc to eventually obtain assocLNorig: SNP  A1  p/OR*  score}
  ss4:=' --bfile '+fnameplinkbin+' --assoc fisher-midp counts --out TMP';
  exitcode:=ExecuteProcess(plinkname,ss4);
  if exitcode<>0 then begin
   writeln('Exit code = ',exitcode,': Run plink --assoc to make TMP* files');
   halt;
  end;
{
write('See TMP* files'); readln;
}
  {Result is TMP.assoc.fisher: CHR  SNP  BP  A1  C_A  C_U  A2  P  OR}
  {Write this to T2.tmp: SNP  A1  p/OR*  score; p/OR = sort variables, score = ln(OR)}
  assign(f1,'TMP.assoc.fisher'); reset(f1);
  readln(f1);	{skip header}
  assign(f2,'T2.tmp');	{SNP  A1  (1-p)/OR*  score=ln(OR)}
  rewrite(f2);
  allsnps:=0; N1:=0;

  {Evaluate all variants}
  while not eof(f1) do begin
   inc(allsnps);
   readstr(f1,sss);	{CHR; not used}
   readstr(f1,sss);	{ 1) SNP }
   write(f2,sss);

   readstr(f1,sss);	{BP; not used}
   readstr(f1,A1);	{ 2) A1 }
   write(f2,bl,A1);

   read(f1,aa,cc);	{integer}
   readstr(f1,sss);	{A2; not used}
   inc(N1,aa+cc);
   bb:=casealleles-aa;
   dd:=ctrlalleles-cc;

   {Write sort variable, either 1-p or OR*; large values are always good}
   read(f1,pp);   {p-value}
   pp:=one-pp;    {for sorting}
   {$I-}
   read(f1,odds); {may be letters in plink output}
   {$I+}
   if IOresult<>0 then odds:=(aa+half)*(dd+half)/( (bb+half)*(cc+half) );
   readln(f1);

   if sortbyOR	  {Sorting to obtain proper subsets of "best" variants}
   then write(f2,bl,odds:8:6)        { 3) sort by OR* }
   else write(f2,bl,pp:14:12);       { 3) sort by 1-p }

   if scorePPV
   then scovar:=(aa+half)/(aa+cc+1)
   else scovar:=(aa+dd)/(aa+bb+cc+dd);

   writeln(f2,bl,scovar:9:6);        { 4) score = PPV; f2=T2.tmp }
  end; {while not eof f1}

  close(f1); close(f2);
{
write('makeLNor: See T2.tmp'); readln;
}
  {Sort T2.tmp by 1-p or OR* into assocLNorig.txt: SNP  A1  p/OR*  score=ln(OR)}
  ss4:=' -o assocLNorig.txt -r -k 3 T2.tmp';
  exitcode:=ExecuteProcess('/usr/bin/sort',ss4);
  if exitcode<>0 then begin
   writeln('Exit code = ',exitcode,': Sort T2.tmp into assocLNorig.txt');
   halt;
  end;
{
write('makeLNor: Check assocLNorig.txt; field 3 should be decreasing'); readln;
}
{
  Now, assocLNorig.txt, length = all SNPs, is the basis for scoring.
  SNP  A1  (1-p)/OR*  score=ln(OR); (1-p) and OR* are decreasing
}
  pp:=N1/( two*allsnps*(pheno[1]+pheno[2]) );	{proportion of A1 alleles}
  writeln(logfile,'Prop. of A1 alleles, P = ',pp:6:4);
  pp:=one-two*pp*(one-pp);
  writeln(logfile,'Random ACC,  1-2P(1-P) = ',pp:6:4);
  writeln(logfile,allsnps,' variants in score file, assocLNorig.txt');

  if paramcount=4 then begin
   writeln(logfile,'Computation stopped due to presence of 4th parameter');
   close(logfile);
   halt;
  end;
 end; {makeLNor}



 PROCEDURE MAKEASSOC;
{
 Do this and later procedures for each set of "best" NNsnp SNPs.
 For given NNsnp best SNPs to use, compute pattern file, assocLN0.txt,
 and score0.profile. Save phenoOrig and scoreOrig for each individual.

 Make initial scores based on NNsnp variants, using assocLN0.txt as
 pattern file; assocLN0 are the first NNsnp lines of assocLNorig.txt.
 --score parameters: filename  SNP_col  A1_col  score_col.
}
 var
  g1,g2:text;
  Nind,kk4:integer;
  ss4:string;

 begin {makeassoc}
  assign(g1,'assocLNorig.txt'); {SNP  A1  p/OR*  score;  p/OR* not used further}
  reset(g1);
  assign(g2,'assocLN0.txt');
  rewrite(g2);

  {Copy NNsnp best SNPs from g1 to g2}
  for kk4:=1 to NNsnp do begin
   readln(g1,ss4); writeln(g2,ss4);
  end;
  close(g1); close(g2);
{
write('subset=',isnp,', makeassoc: Check assocLN0.txt '); readln;
has only NNsnp SNPs, in correct order
}

  {Compute score0 based on NNsnp SNPs => SCORE0.profile; do this maxsnp=13 times}
  ss4:=' --bfile '+fnameplinkbin+' --score assocLN0.txt 1 2 4 --out SCORE0';
  exitcode:=ExecuteProcess(plinkname,ss4);
  if exitcode<>0 then begin
   writeln('Exit code = ',exitcode,': Run plink --score a)');
   halt;
  end;
{
write('subset=',isnp,', makeassoc: Check SCORE0.profile '); readln;
}
{
  Result is SCORE0.profile: FID IID PHENO CNT CNT2 SCORE.
  In SCORE0.profile, mean score for pheno=2 > mean score for pheno=1.
  Keep  FID IID pheno scoreOrig  in memory for all individuals.
}
  assign(g1,'SCORE0.profile'); reset(g1); readln(g1); {skip header}
  Nind:=0;
  while not eof(g1) do begin
   inc(Nind);
   readstr(g1,FID[Nind]); 		{FID}
   readstr(g1,IID[Nind]); 		{IID}
   read(g1,phenoOrig[Nind]);		{pheno}
   for kk4:=1 to 2 do readstr(g1,ss4);	{CNT CNT2; not used}
   readln(g1,scoreOrig[Nind]);		{score = ln(OR)}
  end; {while not eof}
  {writeln('== ',Nind,' pheno and score values read into memory');}
  close(g1); {SCORE0.profile}
  Nindall:=Nind;
 end; {makeassoc}



 procedure analyze1;
{
 Run this for each "leave" individual, whose score is scoreOrig[leave].
 NNsnp = best SNPs to use.
}
 label 44;
 var
  g1,g2:text;
  chrom,ss4,sss:string;
  kk4,Nctrl,caseall,ctrlall:integer;
  mal,phe,predict:shortint; {predict=2 for case; predict=1 for ctrl}
  sscore:vector; {scores for controls}
  rrr,scoremax,prop,tt,ssum,sscoreold,step,aa,bb,cc,dd,pp,odds,scovar1:real;

 begin {analyze1}

  {Run reduced fileset, rem1, to make ASSOC1.assoc.fisher}

  {Make file, remove1.txt, containing "leave" individual to be removed}
  assign(g1,'remove1.txt'); rewrite(g1);
  writeln(g1,FID[leave],bl,IID[leave]);
  close(g1);

  {Make reduced fileset, rem1, without "leave"}
  ss4:=' --bfile '+fnameplinkbin+' --remove remove1.txt --make-bed --out rem1';
  exitcode:=ExecuteProcess(plinkname,ss4);
  if exitcode<>0 then begin
   writeln('Exit code = ',exitcode,': Run plink without individual b) ',leave);
   halt;
  end;

  {Run --assoc to obtain ASSOC1.assoc.fisher, without "leave"}
  {rem1 still contains all SNPs!}
  ss4:=' --bfile rem1 --assoc fisher-midp counts --out ASSOC1';
  exitcode:=ExecuteProcess(plinkname,ss4);
  if exitcode<>0 then begin
   writeln('Exit code = ',exitcode,': Run plink without individual c) ',leave);
   halt;
  end;
  {Result is ASSOC1.assoc.fisher: CHR  SNP  BP  A1  C_A  C_U  A2  P  OR}
{
write('leave=',leave,': Check ASSOC1.assoc.fisher '); readln;
}
{contains all SNPs, in original order, not ordered by p/OR}

  {Read ASSOC1.assoc.fisher; write assocLN1 suitably modified}
  assign(g1,'ASSOC1.assoc.fisher'); reset(g1); readln(g1); {skip header}
  assign(g2,'assocLN1.txt'); rewrite(g2);

  if phenoOrig[leave]=2 then begin
   caseall:=casealleles-2;
   ctrlall:=ctrlalleles;
  end else begin
   caseall:=casealleles;
   ctrlall:=ctrlalleles-2;
  end;

  {Evaluate all variants}
  while not eof(g1) do begin
   readstr(g1,chrom);	{CHR; not used}
   readstr(g1,sss);	{ 1) SNP }
   write(g2,sss);

   readstr(g1,sss);	{BP; not used}
   readstr(g1,sss);	{ 2) A1 }
   write(g2,bl,sss);

   read(g1,aa,cc);	{C_A  C_U}
   readstr(g1,sss);	{A2; not used}
   bb:=caseall-aa;
   dd:=ctrlall-cc;

   {Write sort variable, either 1-p or OR*}
   read(g1,pp);
   pp:=one-pp;
   {$I-}
   read(g1,odds);
   {$I+}
   if IOresult<>0 then odds:=(aa+half)*(dd+half)/( (bb+half)*(cc+half) );
   readln(g1);

   if sortbyOR		{Sorting to obtain proper subsets of "best" variants}
   then write(g2,bl,odds:8:6) 	{ 3) sort by OR* }
   else write(g2,bl,pp:14:12);	{ 3) sort by 1-p }

   if scorePPV
   then scovar1:=(aa+half)/(aa+cc+1)
   else scovar1:=(aa+dd)/(aa+bb+cc+dd);

   writeln(g2,bl,scovar1:9:6);	{ 4) score, PPV }
  end; {while not eof g1}


  close(g1); close(g2);  {assocLN1.txt: SNP  A1  score=ln(OR)}
{
write('leave=',leave,': Check assocLN1.txt, 4 fields '); readln;
}
  {Sort "assocLN1.txt" by p/OR into T1.tmp: SNP  A1  p/OR*  score=ln(OR)}
  ss4:=' -o T1.tmp -r -k 3 assocLN1.txt';
  exitcode:=ExecuteProcess('/usr/bin/sort',ss4);
  if exitcode<>0 then begin
   writeln('Exit code = ',exitcode,': Sort assocLN1.txt into T1.tmp');
   halt;
  end;
{
write('leave=',leave,': Check T1.tmp, all SNPs, 4 fields ordered by decreasing field 3 '); readln;
}
  {Save NNsnp lines of T1.tmp as pattern.txt}
  assign(g1,'T1.tmp'); reset(g1);
  assign(g2,'pattern.txt'); rewrite(g2); {4 fields}
  for kk4:=1 to NNsnp do begin
   readln(g1,ss4); writeln(g2,ss4);
  end; {for kk4}
  close(g1); close(g2);

  {Make scores, TMP.profile, based on the best NNsnp variants}
  {--score filename variantIDcol alleleCol scoreCol}

  ss4:=' --bfile rem1 --score pattern.txt 1 2 4 --out TMP'; {TMP.profile}
 {TMP.profile: FID  IID  PHENO  CNT  CNT2  SCORE, with header line}
  exitcode:=ExecuteProcess(plinkname,ss4);
  if exitcode<>0 then begin
   writeln('Exit code = ',exitcode,': Run plink --score without individual d) ',leave);
   halt;
  end;

  {Copy TMP.profile into TMP2.profile, no header}
  assign(g1,'TMP.profile');  reset(g1); readln(g1);
  assign(g2,'TMP2.profile'); rewrite(g2);
  while not eof(g1) do begin
   readln(g1,ss4); writeln(g2,ss4);
  end; {while}
  close(g1); close(g2);

  {Sort TMP2.profile into SCORE1.profile, sorted by SCORE, increasing}
  ss4:=' -o SCORE1.profile -k 6 TMP2.profile';
  exitcode:=ExecuteProcess('/usr/bin/sort',ss4);
  if exitcode<>0 then begin
   writeln('Exit code = ',exitcode,': Sort TMP.profile into SCORE1.profile');
   halt;
  end;
{
write('Check SCORE1.profile isnp=',iSNP,' leave=',leave,' no header '); readln;
should be increasing
}

  {Find percentile for controls in SCORE1.profile, which should be increasing}

  {Read SCORE1 into memory (sscore) for pheno=1}
  assign(g1,'SCORE1.profile'); reset(g1);
  scoremax:=-9.9E10;
  Nctrl:=0;

  while not eof(g1) do begin
   for kk4:=1 to 2 do readstr(g1,ss4); {FID IID; not used}
   read(g1,phe);
   for kk4:=1 to 2 do readstr(g1,ss4); {CNT CNT2; not used}
   readln(g1,rrr);
   if phe=1 then begin
    inc(Nctrl);
    sscore[Nctrl]:=rrr;
    if rrr>scoremax then scoremax:=rrr;
   end; {if}
  end; {while not eof}

  close(g1);
  if scoremax>999999.0 then begin
   writeln('Some score exceeds 999,999');
   halt;
  end;

  {Sort sscore (contains only controls) by increasing values}
  hpsort(Nctrl,sscore);

  {Find percentile for sscore in controls}
  {mal=1 for 95%, mal=2 for 99%}
  for mal:=1 to ll do begin  {ll = 2 or 3}
   tt:=percent[mal]; {0.95 or 0.99}
   sscoreold:=0.0; ssum:=0.0; step:=one/Nctrl;
   prop:=0.0;
   for kk4:=1 to Nctrl do begin
    ssum:=ssum+sscore[kk4];
    prop:=prop+step;

    if prop >= tt then begin
     rrr := (sscore[kk4]+sscoreold)/two;
     goto 44;
    end; {if}

    sscoreold:=sscore[kk4];
   end; {for k4}

   rrr:=sscoreold;

   {Percentile, rrr, of score for controls; step = pp-ppold}
   44:

   {Save results; "predict" = predicted phenotype}
   if scoreOrig[leave] >= rrr then predict:=2 else predict:=1;
{
                Nsnp      %ile phen pred
    table:array[1..maxsnp,1..2,1..2,1..2] of integer;
            pred=1  pred=2
    ctrl=1    xx      xx
    case=2    xx      xx
}
   inc(table[   isnp,mal,phenoOrig[leave],predict ]);
  end; {for mal}

 end; {analyze1}



 procedure printout;
 var
  j1,j2,j3,j4,a,b,c,d:integer;
  PPV,NPV,ACC:real;
 begin
  {            Nsnp      %ile phen pred}
  {table:array[1..maxsnp,1..3,1..2,1..2] of integer;}
  append(outfile);

  for j1:=1 to maxsnp do begin
   write(outfile,Nsnp[j1]);

   for j2:=1 to ll do begin {j2 corresponds to "mal"; ll = 2 or 3}
    for j3:=2 downto 1 do
    for j4:=2 downto 1 do write(outfile,tab,table[j1,j2,j3,j4]);
    a:=table[j1,j2,2,2];
    b:=table[j1,j2,2,1];
    c:=table[j1,j2,1,2];
    d:=table[j1,j2,1,1];
    if (a+c)=0 then PPV:=-9.9 else PPV:=a/(a+c);
    if (b+d)=0 then NPV:=-9.9 else NPV:=d/(b+d);
    if (a+b+c+d)=0 then ACC:=-9.9 else ACC:=(a+d)/(a+b+c+d);
    write(outfile,tab,PPV:5:3,tab,NPV:5:3,tab,ACC:5:3);
   end; {for j2}

   writeln(outfile);
  end; {for J1}
  close(outfile);
 end; {printout}



 procedure mopup;
{
 Uses the DeleteFile function in sysutils
}
 begin {mopup}
  DeleteFile('ASSOC1.assoc.fisher');
  DeleteFile('ASSOC1.log');
  DeleteFile('assocLN0.txt');
  DeleteFile('assocLN1.txt');
  DeleteFile('assocLNorig.txt');
  DeleteFile('pattern.txt');
  DeleteFile('ppp.bed');
  DeleteFile('ppp.bim');
  DeleteFile('ppp.fam');
  DeleteFile('ppp.log');
  DeleteFile('rem1.bed');
  DeleteFile('rem1.bim');
  DeleteFile('rem1.fam');
  DeleteFile('rem1.log');
  DeleteFile('remove1.txt');
  DeleteFile('SCORE0.log');
  DeleteFile('SCORE0.profile');
  DeleteFile('SCORE1.profile');
  DeleteFile('T1.tmp');
  DeleteFile('T2.tmp');
  DeleteFile('TMP2.profile');
  DeleteFile('TMP.assoc.fisher');
  DeleteFile('TMP.log');
  DeleteFile('TMP.profile');
 end;  {mopup}



begin {L1outPRS}
 {Make ASSOC.assoc.fisher}
 init;
{
 Check values of score
 Write assocLNorig.txt
}
 makeLNor;

 for isnp:=1 to maxsnp do begin
  NNsnp:=Nsnp[isnp]; {best NNsnp variants to use}

  makeassoc;  

  append(reportfile);
  ti(reportfile);
  writeln(reportfile,'  Starting SNP set',isnp:3,' of ',maxsnp,'  Nsnp=',NNsnp);
  close(reportfile);
{
  "leave" = sequential number of the individual left out.
  Analyze Nind-1 individuals
}
  for leave:=1 to Nindall do analyze1;
 end; {for isnp}

 printout;

 writeln;
 writeln(outfilename,' written');
 append(reportfile);
 ti(reportfile);
 writeln(reportfile,'  ',outfilename,' written');
 writeln(reportfile,'All done');
 close(reportfile);

 ti(logfile);
 writeln(logfile,': ',outfilename,' written');
 close(logfile);

 mopup; {remove left-over files}
end.  {L1outPRS}
