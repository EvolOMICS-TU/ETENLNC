#!/usr/bin/perl
package Genecare::PerlOligos;
use Exporter;
use vars qw/@ISA @EXPORT $VERSION/;
use strict;
#user warnings;
@ISA = ('Exporter');
@EXPORT = qw/cal_dimer cal_dimer_2 gc tm/; # subs to export


my ($fpname, $fprimer, $rpname, $rprimer,$d,$pic,$output);

my @result = (); 

my (%oligo_dH, %oligo_dH_full, %oligo_dS, %oligo_dS_full, %genetic_code);

my $oligo_conc = 200; #in nM
my $mg_conc=0.7; #in mM
my $monovalent_cation_conc=50; #in mM
my $dntp_conc=0.2; #in mM

my $repeat = 3;
my $run = 4;
my $numlnc=0;
my $reverse_primer_f;
my $reverse_primer_r;

my $baoding=0;

# Primer-dimer parameters
my $pd_full=0;
my $pd_extensible=1;
my $pd_temperature=37;
my $min_leng=0;
my $sdelta=0;
# More global variables

my $pos_a=0;
my $pos_b=0;

my $number_1=0;
my $number_2=0;
my $number_3=0;
my $number_4=0;
my $number_5=0;
my $number_6=0;

my $score_1=0;
my $score_2=0;

my (
	$primer_f, $primer_r, $pos,
	$rprimer_r,
	$pd, @score_sort, $reverse,
	$pfkeys, $pkeys, @PF, @PR, %primer_hash,
	$gc_exclude, $gc_percent_ex, @primer_pairs,
	$prl, $pfl, $pl, @bind_string, %rating_hash, @score,
	$fprimer_tm, $fprimer_len, $fprimer_ds, $fprimer_dh, $fprimer_dg, $fprimer_gc,
	$rprimer_tm, $rprimer_len, $rprimer_ds, $rprimer_dh, $rprimer_dg, $rprimer_gc,
);

my %oligo_dG=(
	qw(initC 0.98 	initG 0.98 
	initA 1.03 	initT 1.03), 
);
sub cal_dimer {
	($fpname, $fprimer, $rpname, $rprimer,$d,$pic,$output) = @_;# argument list
	load_data();
	recalculate_dG();
	@result=();
	get_tm($fprimer, $rprimer, $pd_full,$pic,$output);
}

sub get_tm {
	my ($report) = @_;
	$fprimer = uc($fprimer);
	$rprimer = uc($rprimer);
		
	my ($deltaG, $deltaH, $deltaS);
	my $oligo_conc_mols = $oligo_conc / 1000000000;
	if ($fprimer && !check_degenerate($fprimer, 1)) {
		($fprimer_tm, $deltaH, $deltaS) = tm($fprimer);
		$fprimer_gc = int(gc($fprimer));
		$fprimer_tm = sprintf("%.2f", $fprimer_tm);
		$deltaG = $deltaH-((273.15+$pd_temperature)*($deltaS/1000));
		$fprimer_ds = sprintf("%.2f", $deltaS);
		$fprimer_dh = sprintf("%.2f", $deltaH);
		$fprimer_dg = sprintf("%.2f", $deltaG);
		$fprimer_len = length($fprimer);
	} else {exit;}	
	if ($rprimer && !check_degenerate($rprimer, 1)) {
		($rprimer_tm, $deltaH, $deltaS) = tm($rprimer);
		$rprimer_gc = int(gc($rprimer));
		$rprimer_tm = sprintf("%.2f", $rprimer_tm);
		$deltaG = $deltaH-((273.15+$pd_temperature)*($deltaS/1000));
		$rprimer_ds = sprintf("%.2f", $deltaS);
		$rprimer_dh = sprintf("%.2f", $deltaH);
		$rprimer_dg = sprintf("%.2f", $deltaG);
		$rprimer_len = length($rprimer);
	}else {exit;}
	print "lncRNA is :$fpname, length : $fprimer_len\t";
	print "mRNA is : $rpname,  length :$rprimer_len\n";
	my $score_1;
  my $score_2;
  my $score_3;
  my $score_4;
  my $score_5;
  my $step;	
	
	if($rprimer_len>=9000)  
  {
    $step= int $rprimer_len/5;
    #print "%%%%%%%%%%%%%step:$step\n";
  	my $zch_1=substr($rprimer,0,$step);
  	primer_dimer($fprimer,$zch_1,1);
    my $pos_1=$rating_hash{$score_sort[0]};
    $score_1=$score_sort[0];

    my $start=$step-$fprimer_len;
    my $end=2*$step;
  	my $zch_2=substr($rprimer,$start,$end);
    primer_dimer($fprimer,$zch_2,1);
    my $pos_2=$rating_hash{$score_sort[0]};
    $score_2=$score_sort[0];

    $start=2*$step-$fprimer_len;
    $end=3*$step;
  	my $zch_3=substr($rprimer,$start,$end);
    primer_dimer($fprimer,$zch_3,1);
    my $pos_3=$rating_hash{$score_sort[0]};
    $score_3=$score_sort[0];

    $start=3*$step-$fprimer_len;
    $end=4*$step;
  	my $zch_4=substr($rprimer,$start,$end);
    primer_dimer($fprimer,$zch_4,1);
    my $pos_4=$rating_hash{$score_sort[0]};
    $score_4=$score_sort[0];

    $start=4*$step-$fprimer_len;
  	my $zch_5=substr($rprimer,$start,$rprimer_len);
    primer_dimer($fprimer,$zch_5,1);
    my $pos_5=$rating_hash{$score_sort[0]};
    $score_5=$score_sort[0];
    if ($score_1<=$score_2 && $score_1<=$score_3 && $score_1<=$score_4&& $score_1<=$score_5)
    {
    	$score_sort[0]=$score_1;  
    	$pos= $pos_1;  
    }elsif($score_2<=$score_1 && $score_2<=$score_3 && $score_2<=$score_4&& $score_2<=$score_5)
    {
    	$score_sort[0]=$score_2; 
    	$pos= $pos_2;   
    }elsif ($score_3<=$score_1 && $score_3<=$score_2 && $score_3<=$score_4 && $score_3<=$score_5)
    {
    	$score_sort[0]=$score_3;
    	$pos= $pos_3;   
    }elsif ($score_4<=$score_1 && $score_4<=$score_3 && $score_4<=$score_2 && $score_4<=$score_5)
    {
    	$score_sort[0]=$score_4;
    	$pos= $pos_4;   
    }else
    {
    	$score_sort[0]=$score_5;
    	$pos= $pos_5;    
    }	
##
#	 if($fprimer_len>=$rprimer_len)
#			{
#				$min_leng=$rprimer_len;
#			}
#			else
#			{
#				$min_leng=$fprimer_len;
#			}
#		  $sdelta=$score_sort[0]/$min_leng;
#		  
#			if ($score_sort[0]<0) {
 #     if($sdelta<=$d)
  #    { 
#			  push (@result,"$fpname\t $fprimer_len\t");
#      	push (@result,"$rpname\t $rprimer_len\t");
#      	push (@result, "$score_sort[0]\t $sdelta\t");
#			  draw_dimer($fprimer, $rprimer, $pos,$score_sort[0],$d);
#        open(PDM,">>$output") or die "Can't open file: $!\n";
#			  print PDM @result;	
#			  close (PDM);
#			}
#			$numlnc++;
#		
#			}
#			return;	
#	}	
	if ($score_sort[0]<0) {
			if($pic eq 'T' || $pic eq 't')
				  {draw_dimer($fprimer, $rprimer, $pos,$score_sort[0],$d);}
				  open(PDM,">>$output") or die "Can't open file: $!\n";
				   print PDM @result;	
				   close (PDM);
				$numlnc++;
				}
				return;	
		}	
		
	my $repeat_real = $repeat-1;
		if ($rprimer && !check_degenerate($rprimer)) {
			
			primer_dimer($fprimer,$rprimer,1);
			$score_1=$score_sort[0];
			$pos_a=$rating_hash{$score_sort[0]};
			$reverse_primer_r=reverse($rprimer);
      primer_dimer($fprimer,$reverse_primer_r,1);
      $score_2=$score_sort[0];
			$pos_b=$rating_hash{$score_sort[0]};
		  if($score_1<=$score_2)
      { 
    	$score_sort[0] = $score_1;
    	$pos=$pos_a;
    	primer_dimer($fprimer,$rprimer,1);
      if ($score_sort[0]<0) {
			draw_dimer($fprimer,$rprimer,$pos,$score_sort[0],$d);
			open(PDM,">>$output") or die "Can't open file: $!\n";
			print PDM @result;	
			close (PDM);
			$numlnc++;
			}
      }
			
		if($score_1>$score_2)
      { 
    	$score_sort[0] = $score_2;
    	$pos=$pos_b;
		if ($score_sort[0]<0) {				
			draw_dimer($fprimer,$reverse_primer_r, $pos,$score_sort[0],$d);
			open(PDM,">>$output") or die "Can't open file: $!";
			print PDM @result;	
			close (PDM);
			$numlnc++;
		
			}
      }		
		}
	
}	
sub cal_dimer_2 {
	($fpname,$fprimer,$rpname,$rprimer,$d,$pic,$output) = @_;# argument list
	load_data();
	recalculate_dG();
	@result = ();
	
	get_tm_2($fpname,$fprimer,$rpname,$rprimer,$d,$pic,$output);
}	
sub get_tm_2 {
	my ($report) = @_;

	$fprimer = uc($fprimer);
	$rprimer = uc($rprimer);

		
	my ($deltaG, $deltaH, $deltaS);
	my $oligo_conc_mols = $oligo_conc / 1000000000;
	if ($fprimer && !check_degenerate($fprimer, 1)) {
		($fprimer_tm, $deltaH, $deltaS) = tm($fprimer);
		$fprimer_gc = int(gc($fprimer));
		$fprimer_tm = sprintf("%.2f", $fprimer_tm);
		$deltaG = $deltaH-((273.15+$pd_temperature)*($deltaS/1000));
		$fprimer_ds = sprintf("%.2f", $deltaS);
		$fprimer_dh = sprintf("%.2f", $deltaH);
		$fprimer_dg = sprintf("%.2f", $deltaG);
		$fprimer_len = length($fprimer);
	}
	else { exit; }
	
	if ($rprimer && !check_degenerate($rprimer, 1)) {
		($rprimer_tm, $deltaH, $deltaS) = tm($rprimer);
		$rprimer_gc = int(gc($rprimer));
		$rprimer_tm = sprintf("%.2f", $rprimer_tm);
		$deltaG = $deltaH-((273.15+$pd_temperature)*($deltaS/1000));
		$rprimer_ds = sprintf("%.2f", $deltaS);
		$rprimer_dh = sprintf("%.2f", $deltaH);
		$rprimer_dg = sprintf("%.2f", $deltaG);
		$rprimer_len = length($rprimer);
	}
	else {exit;}
	
		print "lncRNA is :$fpname, length : $fprimer_len\t";
	  print "mRNA is : $rpname,  length :$rprimer_len\n";
   my $repeat_real = $repeat-1;
		if ($rprimer && !check_degenerate($rprimer)) {
			
			primer_dimer($fprimer,$rprimer,1);
			$score_1=$score_sort[0];
			$pos_a=$rating_hash{$score_sort[0]};
			
			$reverse_primer_r=reverse($rprimer);
      primer_dimer($fprimer,$reverse_primer_r,1);
      $score_2=$score_sort[0];
			$pos_b=$rating_hash{$score_sort[0]};
			

			
		  if($score_1<=$score_2)
      { 
    	$score_sort[0] = $score_1;
    	$pos=$pos_a;
    	primer_dimer($fprimer,$rprimer,1);
     if ($score_sort[0]<0) {
			draw_dimer($fprimer,$rprimer,$pos,$score_sort[0],$d);
			open(PDM,">>$output") or die "Can't open file: $!";
			print PDM @result;	
			close (PDM);
			$numlnc++;
		
		}
      }
			
		if($score_1>$score_2)
      { 
    	$score_sort[0] = $score_2;
    	$pos=$pos_b;
		if ($score_sort[0]<0) {	
		  draw_dimer($fprimer,$reverse_primer_r, $pos,$score_sort[0],$d);
		  open(PDM,">>$output") or die "Can't open file: $!\n";
		  print PDM @result;	
		  close (PDM);
		  $numlnc++;
		}
      }
}
}


sub check_degenerate {
	$_ = shift;
	if (/[^ATGC]/i) {
		print "One of your sequences has a degenerate or non-DNA character,please check\n";
		return 1;
	} 
}


# %GC	
sub gc {
	$_ = $_[0];
	my($gc,$countgc,$counttotal);
	$gc=0;
	$countgc=0;
		
	$countgc = tr/GCgc/GCgc/;
	$counttotal = length();
	
	$gc = $countgc/$counttotal*100;
	return $gc;
}
sub tm {
	my ($primer) = @_;
	$primer = uc($primer); # if user enters primer directly as lower-case
	my ($i, $nn, $initterm, $endterm);
	my $primer_len = length($primer);
	my ($deltaH, $deltaS);
		
	#-----------------------------#
	# calculate deltaH and deltaS #
	#-----------------------------#

	for ($i=0; $i<$primer_len-1; $i++) {
		$nn = substr($primer, $i, 2);
		
		$deltaH+= $oligo_dH{$nn};
		$deltaS+= $oligo_dS{$nn};
	}
		
	#-------------------------#
	# initial term correction #
	#-------------------------#

	$initterm="init" . substr($primer, 0, 1);
	$deltaH+= $oligo_dH{$initterm};
	$deltaS+= $oligo_dS{$initterm};
	
	$endterm="init" . substr($primer, -1, 1);
	$deltaH+= $oligo_dH{$endterm};
	$deltaS+= $oligo_dS{$endterm};
				
	# Tm at 1M NaCl
	# $tm= ($deltaH * 1000) / ($deltaS + (1.987 * log($oligo_conc / 4))) - 273.15;
	
	#------------------------------------------#
	# correct for salt concentration on deltaS #
	#------------------------------------------#
	
	# Big problems if [dNTPs] > [Mg++] !!  This is a quick fix ...
	my $salt_correction;
	if ($mg_conc > $dntp_conc) {
		$salt_correction = sqrt($mg_conc - $dntp_conc);
	} else {
		$salt_correction = 0;
	}
	
	my $na_eq=($monovalent_cation_conc + 120 * $salt_correction)/1000;
	
	# deltaS correction:
	$deltaS += (0.368 * ($primer_len - 1) * log($na_eq));
	
	my $oligo_conc_mols = $oligo_conc / 1000000000;

	# Salt corrected Tm
	# NB - for PCR I'm assuming for the moment that the [strand target] << [oligo]
	# and that therefore the C(t) correction term approx equals [oligo]
	my $corrected_tm=(($deltaH * 1000) / ($deltaS + (1.987 * log($oligo_conc_mols/4)))) - 273.15;
	return ($corrected_tm, $deltaH, $deltaS);
}


sub primer_dimer {
		
	my ($primer_f, $primer_r, $pd_full) = @_;
	return unless ($primer_f) && ($primer_r);
			
	my ($k, $l);
	@score=();
	%primer_hash=();
	@score_sort=();
	@bind_string=();
	%rating_hash=();
		
	# $pl = greatest length
	$pfl=length($primer_f);
	$prl=length($primer_r);
	$pl = ($pfl>$prl ? $pfl : $prl);
	
	my $rcompr = reverse(complement($primer_r));
	my $rcomprlc = lc($rcompr);
	my $fprimer_r=lc(reverse($primer_f));
	$rprimer_r=reverse($primer_r);
	
	# create a binding array for each of the four bases
	for $l (0 .. $pfl-1) {
		my $mbase = substr($fprimer_r, $l, 1);
		$primer_hash{$mbase}[$l]=1;
		for $k (qw/a g c t/) {
			$primer_hash{$k}[$l] ||=0;
		}
	}
		
	# create the primer matrix
	my @primer_comp;
	for $k (0 .. $prl-1) {
		$primer_comp[$k]=$primer_hash{substr($rcomprlc, $k, 1)};
	}
		
	# read each combination from the matrix, calculate dG for each dimer
	my $pd_len = ($pd_full ? $pfl+$prl-1 : $pl-2);
	for $k (0 .. $pd_len) {
		$score[$k]=0;
		my $bind;
		my $score_p=0;
		
		# extensible primer short-circuit - ignore all primers that will
		# not create extensible (i.e. amplifiable) dimers
		my $start = $k>$pfl-1 ? $pfl-1 : $k;
		my $end = $k>$prl-1 ? $prl-1 : $k;
		if ($pd_extensible && !$pd_full) {
			next unless $primer_comp[0][$start] == 1;
			next unless $primer_comp[$end][$start-$k] == 1;
		}
		
		# read the binding data
		for $l (0 .. $prl-1) {
			if (($k-$l)<$pfl) {
				$bind .= $primer_comp[$l][$k-$l] if ($k-$l)>=0;
			} else {
				# spacer
				$bind .= "2";
			}
		}
		
		# Single matched bases surrounded by mismatches are unstable,
		# so we remove them with the regexp (look ahead is needed otherwise
		# strings of consecutive match/mismatches are not caught)
		$bind =~ s/01(?=[^1])/00/gx;
		
		# Short circuit if there's nothing to bind
		next unless $bind =~ /[1]/;
		
		# Find start and end of similarity
		my ($pb_init,$pb_end);
		for $l (0 .. length($bind)-1) {
			# at first I tried finding the initiating terminal bases with
			# regexps, but that was much slower ...
			if (substr($bind, $l, 1) eq "1") {
				defined($pb_init) || ($pb_init = $l);
				$pb_end=$l;
			}
		}
				
		if (defined($pb_init)) {
			# deltaG calculation
			for $l ($pb_init .. $pb_end-1) {
				next if substr($bind, $l, 2) eq "00";
				next if substr($bind, $l, 1) eq "2";
				$score_p+=$oligo_dG{substr($primer_f, $pfl-$k+$l-1, 2).substr($rprimer_r, $l, 2)};
			}
			
			# init term corrections
			my $initterm="init" . substr($rprimer_r, $pb_init, 1);
			$score_p+= $oligo_dG{$initterm};
			
			my $endterm="init" . substr($rprimer_r, $pb_end, 1);
			$score_p+= $oligo_dG{$endterm};
			
			# add to the hash ...
			$score[$k]=sprintf("%.2f",$score_p);
			$bind_string[$k]=$bind;
			$rating_hash{$score[$k]}=$k;
		}
	}
	
	# sort the dimers to give the most stable:	
	@score_sort = sort { $a <=> $b } @score;
		
	# Returns the most stable dimer
	return $score_sort[0];
}


sub draw_dimer {
	# This all seems a bit cumbersome!!
	my ($primer_f, $primer_r, $pos,$score,$d) = @_;
	
	my $rprimer_r=reverse($primer_r);
	my $dimer_binding="";
	my $binding_site="";
	my $binding_site1="";
	my $pr_space="";
	my $fspace="";
	my $rspace="";
	my $tang=0;
	my $shan=0;
	my $tianjin=0;
	my $beijing=0;
	my $hongqiao=0;
	my $nankai=0;
	my $heping=0;
	my $langfang=0;

			
	my $fspace_def = $pl-$pfl>0 ? $pl-$pfl : 0;
	$fspace=" "x($fspace_def+($pos>$pl-1?$pos-$pl+1:0));
	
	if ($pos+1>=$pfl) {
		$rspace=" "x($pl-$pos-1);
	} else {
		$rspace=$fspace;
	}
	
	$pr_space=" "x($pfl-$pos-1);
	##make the space clear
	my $beiyi;
	if(length($fspace)>(length($rspace)+length($pr_space)))
	{
		$beiyi =" "x(length($fspace)-(length($rspace)+length($pr_space)));
	}
	else 
	{
		$beiyi =" "x((length($rspace)+length($pr_space))-length($fspace));
	}

	##end
	for my $j (0 .. $pos) {
		next unless $j < $prl;
		if (substr($bind_string[$pos],$j,1)==1) {
			$dimer_binding=$dimer_binding."|"
		} elsif (substr($bind_string[$pos],$j,1)==0) {
			$dimer_binding=$dimer_binding."."
		} else {
			$dimer_binding=$dimer_binding." "
		}
	}
 $tang=$dimer_binding =~ tr/ / /;  
 $shan=$pr_space =~ tr/ / /;
 $tianjin=$rspace =~ tr/ / /;
 $beijing=$tang+$shan+$tianjin;  
 $langfang=$tianjin+$shan;
 
 
 
 $number_1=length($fspace);
 $number_2=$beijing;
 $number_3=$number_2-$number_1+1;
 $hongqiao=length("$fspace"."5' "."$primer_f");
 $nankai=length("$rspace"."   "."$pr_space"."$dimer_binding");
 $number_4=length($primer_f)-($hongqiao-$nankai);
 
 $number_5=$beijing-$langfang+1;
 $heping=length("$rspace"."$pr_space"."3' "."$rprimer_r");
 $nankai=length("$rspace"."   "."$pr_space"."$dimer_binding");
 $number_6=length($rprimer_r)-($heping-$nankai);
 
## push (@result, "$number_3\t");
## push (@result, "$number_4\t");
## push (@result, "$number_5\t");
## push (@result, "$number_6\n");
 
 ##20151026 nDG = score/min_leng
my @data;
@data = split /\|/, $dimer_binding;
shift @data; pop @data; 
$binding_site1 = join '|', @data;
$binding_site = "|".$binding_site1."|";
$sdelta=$score/length($binding_site);
$sdelta = sprintf("%.4f",$sdelta);

if(($pic eq 'T' || $pic eq 't')&&($sdelta<$d))
	{ 
    push (@result,"$fpname\t $fprimer_len\t");
    push (@result,"$rpname\t $rprimer_len\t");
    push (@result, "$score\t $sdelta\t");
		push (@result, "$number_3\t");
		push (@result, "$number_4\t");
		push (@result, "$number_5\t");
		push (@result, "$number_6\n");
		if(length($fspace)>=(length($rspace)+length($pr_space)))
		{
			push (@result, "$beiyi"."5' "."$primer_f"." 3'\n".
		"   "."$dimer_binding\n".
		"3' "."$rprimer_r"." 5'\n\n");
		}
		if(length($fspace)<(length($rspace)+length($pr_space)))
		{
			push (@result, "5' "."$primer_f"." 3'\n".
		"$beiyi"."   "."$dimer_binding\n".
		"$beiyi"."3' "."$rprimer_r"." 5'\n\n");
		}
	}
if(($pic eq 'F' || $pic eq 'f')&&($sdelta<$d))
	{ 
		push (@result,"$fpname\t $fprimer_len\t");
        push (@result,"$rpname\t $rprimer_len\t");
        push (@result, "$score\t $sdelta\t");
		push (@result, "$number_3\t");
		push (@result, "$number_4\t");
		push (@result, "$number_5\t");
		push (@result, "$number_6\n");
	} 
}	

sub complement {
	$_ = shift;
	tr/AGCTagct/TCGAtcga/;
	return $_;
}
sub recalculate_dG {
	# because dG = dH - TdS, and dS is dependent on the salt concentration ...
	my $temperature = shift || $pd_temperature;
	
	# Big problems if [dNTPs] > [Mg++] !!  This is a quick fix ...
	my $salt_correction;
	if ($mg_conc > $dntp_conc) {
		$salt_correction = sqrt($mg_conc - $dntp_conc);
	} else {
		$salt_correction = 0;
	}
	
	my $na_eq=($monovalent_cation_conc + 120 * $salt_correction)/1000;
	
	# the length of each NN dimer is 2, therefore the modifier is 1
	my $entropy_adjust = (0.368 * log($na_eq));
		
	foreach my $key (keys(%oligo_dH_full)) {
		next if $key =~ /init/; # the length of each monomer is 1, thus the modifier of dS is 0 and the values are precalulated
		
		my $dS = $oligo_dS_full{$key} + $entropy_adjust;
		my $dG = $oligo_dH_full{$key}-((273.15+$temperature)*($dS/1000));
		$oligo_dG{$key} = $dG;
	}
}
sub load_data {	
	%oligo_dH=qw(
		AA -7.9 TT -7.9 
		AT -7.2 TA -7.2 
		CA -8.5 TG -8.5 
		GT -8.4 AC -8.4 
		CT -7.8 AG -7.8 
		GA -8.2 TC -8.2 
		CG -10.6 GC -9.8 
		GG -8.0 CC -8.0 
		initC 0.1 initG 0.1 
		initA 2.3 initT 2.3
	);
	
	%oligo_dH_full=(
		qw(AATT -7.9 	TTAA -7.9 
		ATTA -7.2 	TAAT -7.2 
		CAGT -8.5 	TGAC -8.5 
		GTCA -8.4 	ACTG -8.4 
		CTGA -7.8 	AGTC -7.8 
		GACT -8.2 	TCAG -8.2 
		CGGC -10.6 	GCCG -9.8 
		GGCC -8.0 	CCGG -8.0
			
		initC 0.1 	initG 0.1 
		initA 2.3 	initT 2.3),
		
		# Like pair mismatches 
			
		qw(AATA 1.2 	ATAA 1.2
		CAGA -0.9 	AGAC -0.9
		GACA -2.9 	ACAG -2.9
		TAAA 4.7 	AAAT 4.7 
		
		ACTC 0.0 	CTCA 0.0 
		CCGC -1.5 	CGCC -1.5
		GCCC 3.6 	CCCG 3.6 
		TCAC 6.1 	CACT 6.1 
		
		AGTG -3.1 	GTGA -3.1
		CGGG -4.9 	GGGC -4.9
		GGCG -6.0 	GCGG -6.0
		TGAG 1.6 	GAGT 1.6 
		
		ATTT -2.7 	TTTA -2.7
		CTGT -5.0 	TGTC -5.0
		GTCT -2.2 	TCTG -2.2
		TTAT 0.2 	TATT 0.2  ),
		
		# G.T mismatches 
		
		qw(AGTT 1.0  	TTGA 1.0
		ATTG  -2.5 	GTTA  -2.5
		CGGT  -4.1 	TGGC  -4.1
		CTGG  -2.8 	GGTC  -2.8
		GGCT  3.3 	TCGG  3.3
		GGTT  5.8 	TTGG  5.8
		GTCG  -4.4 	GCTG  -4.4
		GTTG  4.1 	GTTG  4.1
		TGAT  -0.1 	TAGT  -0.1
		TGGT  -1.4 	TGGT  -1.4
		TTAG  -1.3 	GATT  -1.3), 
		
		# G.A mismatches 
		
		qw(AATG  -0.6 	GTAA  -0.6
		AGTA  -0.7 	ATGA  -0.7
		CAGG  -0.7 	GGAC  -0.7
		CGGA  -4.0 	AGGC  -4.0
		GACG  -0.6 	GCAG  -0.6
		GGCA  0.5 	ACGG  0.5
		TAAG  0.7 	GAAT  0.7
		TGAA  3.0 	AAGT  3.0), 
		
		# C.T mismatches 
		
		qw(ACTT  0.7 	TTCA  0.7
		ATTC  -1.2 	CTTA  -1.2
		CCGT  -0.8 	TGCC  -0.8
		CTGC  -1.5 	CGTC  -1.5
		GCCT  2.3 	TCCG  2.3 
		GTCC  5.2 	CCTG  5.2 
		TCAT  1.2 	TACT  1.2 
		TTAC  1.0 	CATT  1.0), 
		
		# A.C mismatches 
		
		qw(AATC  2.3	CTAA  2.3
		ACTA  5.3 	ATCA  5.3 
		CAGC  1.9 	CGAC  1.9 
		CCGA  0.6 	AGCC  0.6 
		GACC  5.2 	CCAG  5.2 
		GCCA  -0.7 	ACCG  -0.7
		TAAC  3.4  	CAAT  3.4 
		TCAA  7.6 	AACT  7.6),
	
	);
	
	#--------------------#
	# deltaS (cal/K.mol) #
	#--------------------#
	
	%oligo_dS=qw(
		AA -22.2 TT -22.2 
		AT -20.4 TA -21.3 
		CA -22.7 TG -22.7 
		GT -22.4 AC -22.4 
		CT -21.0 AG -21.0 
		GA -22.2 TC -22.2 
		CG -27.2 GC -24.4 
		GG -19.9 CC -19.9 
		initC -2.8 initG -2.8 
		initA 4.1 initT 4.1 
		sym -1.4
	);
	
	%oligo_dS_full=(
		qw(AATT -22.2 	TTAA -22.2 
		ATTA -20.4 	TAAT -21.3 
		CAGT -22.7 	TGAC -22.7 
		GTCA -22.4 	ACTG -22.4 
		CTGA -21.0 	AGTC -21.0 
		GACT -22.2 	TCAG -22.2 
		CGGC -27.2 	GCCG -24.4 
		GGCC -19.9 	CCGG -19.9
			
		initC -2.8 	initG -2.8 
		initA 4.1 	initT 4.1
		sym -1.4),
		
		# Like pair mismatches
			
		qw(AATA 1.7 	ATAA 1.7
		CAGA -4.2 	AGAC -4.2 
		GACA -9.8 	ACAG -9.8 
		TAAA 12.9 	AAAT 12.9 
		
		ACTC -4.4 	CTCA -4.4 
		CCGC -7.2 	CGCC -7.2 
		GCCC 8.9 	CCCG 8.9 
		TCAC 16.4 	CACT 16.4 
		
		AGTG -9.5 	GTGA -9.5 
		CGGG -15.3 	GGGC -15.3
		GGCG -15.8 	GCGG -15.8
		TGAG 3.6 	GAGT 3.6 
		
		ATTT -10.8 	TTTA -10.8
		CTGT -15.8 	TGTC -15.8
		GTCT -8.4 	TCTG -8.4 
		TTAT -1.5 	TATT -1.5),
		
		# G.T mismatches
		
		qw(AGTT 0.9 	TTGA 0.9
		ATTG  -8.3 	GTTA  -8.3
		CGGT  -11.7 	TGGC  -11.7
		CTGG  -8.0 	GGTC  -8.0
		GGCT  10.4 	TCGG  10.4
		GGTT  16.3 	TTGG  16.3
		GTCG  -12.3 	GCTG  -12.3
		GTTG  9.5 	GTTG  9.5
		TGAT  -1.7 	TAGT  -1.7
		TGGT  -6.2 	TGGT  -6.2
		TTAG  -5.3 	GATT  -5.3), 
		
		# G.A mismatches
		
		qw(AATG  -2.3 	GTAA  -2.3
		AGTA  -2.3 	ATGA  -2.3
		CAGG  -2.3 	GGAC  -2.3
		CGGA  -13.2 	AGGC  -13.2
		GACG  -1.0 	GCAG  -1.0
		GGCA  3.2 	ACGG  3.2
		TAAG  0.7 	GAAT  0.7
		TGAA  7.4 	AAGT  7.4), 
		
		# C.T mismatches
		
		qw(ACTT  0.2 	TTCA  0.2
		ATTC  -6.2 	CTTA  -6.2
		CCGT  -4.5 	TGCC  -4.5
		CTGC  -6.1 	CGTC  -6.1
		GCCT  5.4 	TCCG  5.4 
		GTCC  13.5 	CCTG  13.5
		TCAT  0.7 	TACT  0.7 
		TTAC  0.7 	CATT  0.7), 
		
		# A.C mismatches
		
		qw(AATC  4.6 	CTAA  4.6
		ACTA  14.6 	ATCA  14.6
		CAGC  3.7 	CGAC  3.7 
		CCGA  -0.6 	AGCC  -0.6
		GACC  14.2 	CCAG  14.2
		GCCA  -3.8 	ACCG  -3.8
		TAAC  8.0  	CAAT  8.0 
		TCAA  20.2 	AACT  20.2),
	
	);
	
	
	# Genetic code hash
	%genetic_code=qw(
			TTT F TTC F TTA L TTG L
			CTT L CTC L CTA L CTG L
			ATT I ATC I ATA I ATG M
			GTT V GTC V GTA V GTG V
			TCT S TCC S TCA S TCG S
			CCT P CCC P CCA P CCG P
			ACT T ACC T ACA T ACG T
			GCT A GCC A GCA A GCG A
			TAT Y TAC Y TAA * TAG *
			CAT H CAC H CAA Q CAG Q
			AAT N AAC N AAA K AAG K
			GAT D GAC D GAA E GAG E
			TGT C TGC C TGA * TGG W
			CGT R CGC R CGA R CGG R
			AGT S AGC S AGA R AGG R
			GGT G GGC G GGA G GGG G
	);
}

