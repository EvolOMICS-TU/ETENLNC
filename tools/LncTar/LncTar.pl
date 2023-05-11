#!/usr/bin/perl -w 
use Genecare::PerlOligos;
use strict;
use warnings;
use Getopt::Std;
use List::Util qw[min max];

my %opt=();
my $lncprimerfile='';
my $mRNAprimerfile='';
my $inputfile='';
my $delta=0;
my $kind=0;
my @res = ();

my $pdm_result = "";
my $whetherpic='';
my ($fprimer,$rprimer);
my @primers;
my @temp;
my @temp1;
my @para;
my @pdm_result;
my $list_len = 0;
my $lncPrimerCount=0;

my $usage_1="Usage: perl LncTar.pl [-p kind] [-l lncfile] [-m mrnafile] [-d delta] [-s WhetherPic ] [-o outfile] \n";
my $usage_2="Usage: perl LncTar.pl [-p kind] [-f filename] [-d delta] [-s WhetherPic]  [-o outfile] \n";

if(!(@ARGV==12 || @ARGV==10))
{
   print "param error ,please check out and run again\n";
   print "$usage_1";
   print "$usage_2";
   exit;
}


if(@ARGV==10)
{ 
   if($ARGV[0] eq '-p')
    { $para[1]=$ARGV[1];}
   if($ARGV[2] eq '-p')
    { $para[1]=$ARGV[3];}
   if($ARGV[4] eq '-p')
    { $para[1]=$ARGV[5];}
   if($ARGV[6] eq '-p')
    { $para[1]=$ARGV[7];}
   if($ARGV[8] eq '-p')
    { $para[1]=$ARGV[9];}   
    if($ARGV[0] eq '-f')
    { $para[3]=$ARGV[1];}
    if($ARGV[2] eq '-f')
    { $para[3]=$ARGV[3];}
    if($ARGV[4] eq '-f')
    { $para[3]=$ARGV[5];}
    if($ARGV[6] eq '-f')
    { $para[3]=$ARGV[7];}
    if($ARGV[8] eq '-f')
    { $para[3]=$ARGV[9];}   
    if($ARGV[0] eq '-s')
    { $para[5]=$ARGV[1];}
    if($ARGV[2] eq '-s')
    { $para[5]=$ARGV[3];}
    if($ARGV[4] eq '-s')
    { $para[5]=$ARGV[5];}
    if($ARGV[6] eq '-s')
    { $para[5]=$ARGV[7];}  
    if($ARGV[8] eq '-s')
    { $para[5]=$ARGV[9];}      
    if($ARGV[0] eq '-d')
    { $para[7]=$ARGV[1];}
    if($ARGV[2] eq '-d')
    { $para[7]=$ARGV[3];}
    if($ARGV[4] eq '-d')
    { $para[7]=$ARGV[5];}
    if($ARGV[6] eq '-d')
    { $para[7]=$ARGV[7];}
    if($ARGV[8] eq '-d')
    { $para[7]=$ARGV[9];}   
    if($ARGV[0] eq '-o')
    { $para[9]=$ARGV[1];}
    if($ARGV[2] eq '-o')
    { $para[9]=$ARGV[3];}
    if($ARGV[4] eq '-o')
    { $para[9]=$ARGV[5];}
    if($ARGV[6] eq '-o')
    { $para[9]=$ARGV[7];}
    if($ARGV[8] eq '-o')
    { $para[9]=$ARGV[9];}
}
if(@ARGV==12)
{ 
   if($ARGV[0] eq '-p')
    { $para[1]=$ARGV[1];}
   if($ARGV[2] eq '-p')
    { $para[1]=$ARGV[3];}
   if($ARGV[4] eq '-p')
    { $para[1]=$ARGV[5];}
   if($ARGV[6] eq '-p')
    { $para[1]=$ARGV[7];}
   if($ARGV[8] eq '-p')
    { $para[1]=$ARGV[9];}
   if($ARGV[10] eq '-p')
    { $para[1]=$ARGV[11];}  
    if($ARGV[0] eq '-l')
    { $para[3]=$ARGV[1];}
    if($ARGV[2] eq '-l')
    { $para[3]=$ARGV[3];}
    if($ARGV[4] eq '-l')
    { $para[3]=$ARGV[5];}
    if($ARGV[6] eq '-l')
    { $para[3]=$ARGV[7];}
    if($ARGV[8] eq '-l')
    { $para[3]=$ARGV[9];}
   if($ARGV[10] eq '-l')
    { $para[3]=$ARGV[11];} 
    if($ARGV[0] eq '-m')
    { $para[5]=$ARGV[1];}
    if($ARGV[2] eq '-m')
    { $para[5]=$ARGV[3];}
    if($ARGV[4] eq '-m')
    { $para[5]=$ARGV[5];}
    if($ARGV[6] eq '-m')
    { $para[5]=$ARGV[7];}
    if($ARGV[8] eq '-m')
    { $para[5]=$ARGV[9];}
   if($ARGV[10] eq '-m')
    { $para[5]=$ARGV[11];}  
    if($ARGV[0] eq '-d')
    { $para[7]=$ARGV[1];}
    if($ARGV[2] eq '-d')
    { $para[7]=$ARGV[3];}
    if($ARGV[4] eq '-d')
    { $para[7]=$ARGV[5];}
    if($ARGV[6] eq '-d')
    { $para[7]=$ARGV[7];}
    if($ARGV[8] eq '-d')
    { $para[7]=$ARGV[9];}
   if($ARGV[10] eq '-d')
    { $para[7]=$ARGV[11];}
     if($ARGV[0] eq '-s')
    { $para[9]=$ARGV[1];}
    if($ARGV[2] eq '-s')
    { $para[9]=$ARGV[3];}
    if($ARGV[4] eq '-s')
    { $para[9]=$ARGV[5];}
    if($ARGV[6] eq '-s')
    { $para[9]=$ARGV[7];}
    if($ARGV[8] eq '-s')
    { $para[9]=$ARGV[9];}
   if($ARGV[10] eq '-s')
    { $para[9]=$ARGV[11];}
     if($ARGV[0] eq '-o')
    { $para[11]=$ARGV[1];}
    if($ARGV[2] eq '-o')
    { $para[11]=$ARGV[3];}
    if($ARGV[4] eq '-o')
    { $para[11]=$ARGV[5];}
    if($ARGV[6] eq '-o')
    { $para[11]=$ARGV[7];}
    if($ARGV[8] eq '-o')
    { $para[11]=$ARGV[9];}
   if($ARGV[10] eq '-o')
    { $para[11]=$ARGV[11];}
}
if(!($para[1]==1 ||$para[1]==2))
{
	print "the first param should be 1 or 2!\n";
	exit;
}
if($para[1]==1)    
{
	if(@ARGV!=12)
	{
		 print "param error.\n";
		 print "$usage_1";
		 exit;
  }
        $lncprimerfile=$para[3];
        $mRNAprimerfile=$para[5];
        $delta=$para[7];
        $whetherpic=$para[9];
        $pdm_result=$para[11];
   if((!$lncprimerfile)||(!$mRNAprimerfile)||(!$delta) ||(!$whetherpic)||(!$pdm_result))
   { print "param error,please check carefully.\n";
		 print "$usage_1";
		 exit;
  }
        
}
if($para[1]==2)  
{
	if(@ARGV!=10)
	{
		 print "param error.\n";
		 print "$usage_2";
		 exit;
        }
        $lncprimerfile=$para[3];               
        $whetherpic=$para[5];
        $delta=$para[7];
        $pdm_result=$para[9];
        if((!$lncprimerfile)||(!$delta) ||(!$whetherpic)||(!$pdm_result))
   { print "param error,please check carefully.\n";
		 print "$usage_2";
		 exit;
  }
        
}
open (LOGF, ">$pdm_result")||die"can not open pdm_result!\n"; 
close LOGF;     #clear the file of pdm_result

if($para[1]==1)
{
	      push (@res,"Query\t\t Length_Query\t\t Target\t\t Length_Target\t\t dG\t\t ndG\t\t Start_Position_Query\t\t End_Position_Query\t\t Start_Position_Target\t\t End_Position_Target\n");
	      open(PDM,">>$pdm_result") or die "Can't open file $pdm_result : No such file or directory!\n";
			  print PDM @res;	
			  close (PDM); 
#open (PRIMER_FILE, "<$lncprimerfile") or die "Can't open file $lncprimerfile : No such file or directory!\n";
##
##
open(LNCRNA,"<$lncprimerfile")||die "Can't open $lncprimerfile";
open(LNCTEMP,">lnc_temp.txt")||die "nonononono";
my $flag = 1;
my @temp_lnc;
my @temp_m;
while(<LNCRNA>)
{
	my $temp_line = $_;
	chomp($temp_line);
	$temp_line =~ s/(^\s+|\s+$)//g;#È¥³ýÁ½¶Ë¿Õ°×
	$temp_line =~ s/[\r]$//g;
	if ($temp_line=~/^>/)
	{
		@temp_lnc =split(/\s+/,$temp_line);
		if($temp_lnc[1])
		{
			print "the format of file is error, Please check! The name and the sequence do not put in a row!\n ";
			exit;
		}
		if($flag==1)
		{
			$temp_line = $temp_line."\n";
			$flag=$flag+1;
		}
		else
		{
			$temp_line = "\n".$temp_line."\n";
		}
	}
	elsif ($flag==1)
	{
		print "the format of file is error,maybe you didn't input the name of the sequence or the name didn't start with \">\"!\n ";
		exit;
	}
	print LNCTEMP $temp_line;
}
close LNCRNA;
close LNCTEMP;

open(MRNA,"<$mRNAprimerfile")||die "Can't open $mRNAprimerfile";
open(MRNATEMP,">mRNA_temp.txt")||die "nonononono";
my $flag1 = 1;
while(<MRNA>)
{
	my $temp_line1 = $_;
	chomp($temp_line1);
	$temp_line1 =~ s/(^\s+|\s+$)//g;
	$temp_line1 =~ s/[\r]$//g;
	if ($temp_line1=~/^>/)
	{
		@temp_m =split(/\s+/,$temp_line1);
		if($temp_m[1])
		{
			print "the format of file is error,Please check! The name and the sequence do not put in a row!\n ";
			exit;
		}
		if($flag1==1)
		{
			$temp_line1 = $temp_line1."\n";
			$flag1=$flag1+1;
		}
		else
		{
			$temp_line1 = "\n".$temp_line1."\n";
		}
	}
	elsif ($flag1==1)
	{
		print "the format of file is error,maybe you didn't input the name of the sequence!\n ";
		exit;
	}
	print MRNATEMP $temp_line1;
}
close MRNA;
close MRNATEMP;

##
##
open (PRIMER_FILE, "<lnc_temp.txt") or die "Can't open file $lncprimerfile : No such file or directory!\n";
while (<PRIMER_FILE>) {	   
		$list_len=0;
		my $line = $_;
		chomp($line);
		#$line =~ s/\s+$//;
		$line =~ s/(^s+|s+$)//g;#clear the space both
		if($line =~ /^>/)
		{
			$temp[0]=substr($line,1);#save the name of lncRNA
		}
		next if($line =~ /^>/);
		if($line =~ /^(A|T|C|G|a|t|c|g)/){
			$temp[1]=$line;
			$temp[1]=~s///g; 
			$temp[1]=~s/ //g; 
			$temp[1]=~s/u/t/g;
			$temp[1]=~s/U/T/g;
			$temp[1]=~s/\n//g; 
			push (@primers, [@temp]); 
			$list_len++;
			open (MRNA_FILE,"<mRNA_temp.txt") or die "Can't open file $mRNAprimerfile : No such file or directory!\n";
			while(<MRNA_FILE>)
			{
				my $line1 = $_;
				chomp($line1);
				$line1 =~ s/\s+$//;
				if($line1 =~ /^>/)
				{
					$temp1[0]=substr($line1,1);
				}
				next if($line1 =~ /^>/);
				if($line1 =~ /^(A|T|C|G|a|t|c|g)/){
					$temp1[1]=$line1;
					$temp1[1]=~s///g;
					$temp1[1]=~s/u/t/g;
					$temp1[1]=~s/U/T/g;
					$temp1[1]=~s/ //g; 
					$temp1[1]=~s/\n//g;
					push (@primers, [@temp1]); 
					$list_len++;
				}
				else{
					print "the format of file (-m) is error,Please check.\n ";
					exit;
				}
			}
}

	my $index=$lncPrimerCount/50;
	my $nowCount=0;
	while($index>0)
	{
		$index--;
		$nowCount++;
	}
	if($index<0)
	{
		$nowCount--;
	}
	for (my $i = 0; $i <=0; $i++) {   
			for (my $j =$i+1; $j <= $list_len-1; $j++) {
				cal_dimer($primers[$i][0], $primers[$i][1], $primers[$j][0], $primers[$j][1],$delta,$whetherpic,$pdm_result);
				}	
		}
	$lncPrimerCount++;
	@temp=();
	@primers=();
	close MRNA_FILE;	
  }
close PRIMER_FILE;
}  




if($para[1]==2)
{
        push (@res,"Query\t\t Length_Query\t\t Target\t\t Length_Target\t\t dG\t\t ndG\t\t Start_Position_Query\t\t End_Position_Query\t\t Start_Position_Target\t\t End_Position_Target\n");
	      open(PDM,">>$pdm_result") or die "Can't open file: $!";
			  print PDM @res;	
			  close (PDM);
open(PRIMER_FILE, "<$lncprimerfile") or die "Can't open file $lncprimerfile : No such file or directory!\n";;
while (<PRIMER_FILE>) {	
	$list_len=0;
        my $line = $_;
	      my @temp = split (/[\t ]/, $line);
	      	if($temp[1] eq '' ||$temp[1] eq '' )
  {
  	print "the format of file (-f) is error,Please check.\n ";
  	exit;
  }  
        $temp[1]=~s///g; 
        $temp[1]=~s/ //g; 
        $temp[1]=~s/u/t/g;
		    $temp[1]=~s/U/T/g;
        $temp[3]=~s///g;
        $temp[3]=~s/ //g;
        $temp[3]=~s/u/t/g;
		    $temp[3]=~s/U/T/g;
        $temp[1]=~s/\n//g;
        $temp[3]=~s/\n//g;                       
	      
	cal_dimer_2($temp[0],$temp[1],$temp[2],$temp[3],$delta,$whetherpic,$pdm_result);		
	$lncPrimerCount++;
	@temp=();
	@primers=();

  }
close PRIMER_FILE;
}
close LOGF;
