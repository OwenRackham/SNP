#!/usr/bin/perl -w

use strict;
use warnings;
use DBI;
use lib '/home/rackham/modules';
use rackham;
use Data::Dumper;
use URI::Escape;
use Math::CDF qw(:all);
use Statistics::Basic qw(:all);
use List::Util qw[min max];


my $pheno = $ARGV[0];

		my	$dbh = rackham::DBConnect('superfamily');
	
open FILE, "/home/ephas/Group/KOREAN/SF.csv" or die $!;
print "file opened\n";
my %Ms;#{per person}{per PO} -> how many domains with fathmm score < 1
my %Ns;#{per person} -> how many domains
my %Ks;#{per person}{per PO} -> how many domains
my %Xs;#{per person}{per PO} -> how many with fathmm score < 1
my $j=0;
while (<FILE>) {
	    
	    
		unless ($_ =~ /^Exome/){
		my @line = split(/\t/,$_);
		my (  $sth );
		#print "Before query\n";
		$sth =   $dbh->prepare( "select distinct model.sf,PO_mapping.po,all_score from PO_mapping,model where model.sf = PO_mapping.id and PO_mapping.inherited_from !='' and model.model = $line[3] and PO_mapping.obo in ('$pheno');" );
        $sth->execute;
		#print "After query\n";
		my $i=0;
        while (my ($sf,$PO,$all_score)= $sth->fetchrow_array ) {
			#print "$line[0]\t$model\t$PO\t$all_score\t$line[4]\t$sf\n";
			if($line[4]<1){
				if(exists($Ms{$line[0]}{$sf})){
					$Ms{$line[0]}{$sf} = $Ms{$line[0]}{$sf} + 1;
				}else{
					$Ms{$line[0]}{$sf} = 1;
				}
				if(exists($Xs{$line[0]}{$PO}{$sf})){
					$Xs{$line[0]}{$PO}{$sf} = $Xs{$line[0]}{$PO}{$sf} + 1;
				}else{
					$Xs{$line[0]}{$PO}{$sf} = 1;
				}
			}
			if(exists($Ns{$line[0]}{$sf})){
					$Ns{$line[0]}{$sf} = $Ns{$line[0]}{$sf} + 1;
				}else{
					$Ns{$line[0]}{$sf} = 1;
				}
			if(exists($Ks{$line[0]}{$PO}{$sf})){
					$Ks{$line[0]}{$PO}{$sf} = $Ks{$line[0]}{$PO}{$sf} + 1;
				}else{
					$Ks{$line[0]}{$PO}{$sf} = 1;
				}	
			
			print "$i\n" if(!(++$i%100));
		}		
		print "$line[3]:\t$i\n";
	}
}

my %score;#{per person}{per PO} -> binomial score
my %y;#{per person}{per PO} -> pvalue
my %q;#{per person}{per PO} -> fdr
foreach my $person (keys %Xs){
	my $ns=scalar(keys %{$Ns{$person}});
	foreach my $po (keys %{$Xs{$person}}){
		my $ms=scalar(keys %{$Ms{$person}});
		my $ks=scalar(keys %{$Ks{$person}{$po}});
		my $xs=scalar(keys %{$Xs{$person}{$po}});
		my ($score,$y)=hypergeo_enrich($xs,$ks,$ms,$ns);
		$score{$person}{$po}=$score;
		$y{$person}{$po}=$y;
	}
	my %q_tmp=ConvertP2FDR_1D(\%{$y{$person}});
	$q{$person}=\%q_tmp;
}

my $sth =   $dbh->prepare( "select name from PO_info where obo='$pheno' and po=?");


foreach my $person (keys %y){
	open(PERSON,">$person.$pheno");
	my $ns=scalar(keys %{$Ns{$person}});
	foreach my $po (sort {$score{$person}{$b} <=> $score{$person}{$a}} keys %{$y{$person}}){
		my $ms=scalar(keys %{$Ms{$person}});
		my $ks=scalar(keys %{$Ks{$person}{$po}});
		my $xs=scalar(keys %{$Xs{$person}{$po}});
		if($xs>=3){
			$sth->execute($po);
			my $name=$sth->fetchrow_array;
			print PERSON "$person\t$po\t$name\t$xs\t$ms\t$ks\t$ns\t$score{$person}{$po}\t$y{$person}{$po}\t$q{$person}{$po}\n";
		}
	}
}


sub hypergeo_enrich {
   my $x   = shift;
   my $k   = shift;
   my $m   = shift;
   my $n   = shift;

   my $lscore=0;
   my $y=1;
	
	print "$x\t$k\t$m\t$n\n";
	
   my $p=$k/$n;
   my $exp=$m*$p; # Expected
   my $var=$m*$p*(1-$p); # Variances
   if($var!=0){
       $lscore=($x-$exp)/sqrt($var); # binomial score
       $y=1-pbinom($x, $m, $p); # probability corresponding to fisher's exact test
   }

   $lscore=sprintf "%.2f", $lscore;
   return ($lscore,$y);
}

sub ConvertP2FDR_1D {
       my %y = %{(shift)}; # hash: key (ij), value (p-value)

       my %fdr; # hash: key (ij), value (fdr)

       my (@tmpi,@tmpk);
       foreach my $supra_ij (keys %y){
           push @tmpi,$supra_ij;
           push @tmpk,$y{$supra_ij};
       }

   print "\t\tsort the array by value in an ascending manner\n";
   # sort the array by value in an ascending manner
   my @index_order=sort { $tmpk[$a] <=> $tmpk[$b] } 0 .. $#tmpk;

   print "\t\tcalculate fdr using step-up procedures controlling FDR (BH, 1995)\n";
   # calculate fdr using step-up procedures controlling FDR (BH, 1995)
   my @loc;
   for my $i ( 0 .. $#index_order ) {
       $loc[$i]=min(($#index_order+1)/($i+1)*$tmpk[$index_order[$i]],1);
   }

   my $line=0;
   my $numLines=$#index_order+1; # total number of lines
   my $min_current=1;
   for(my $i=$#index_order;$i>=0;$i--){

       if($loc[$i]<$min_current){
           $min_current=$loc[$i];
       }else{
           $loc[$i]=$min_current;
       }

       # 100 percentage in progress
       #progress_indication($numLines,++$line);

   }

   print "\t\tconvert back into the original order\n";
   for my $i ( 0 .. $#index_order ) {
       $fdr{$tmpi[$index_order[$i]]}=$loc[$i];
       }

       return %fdr;
}