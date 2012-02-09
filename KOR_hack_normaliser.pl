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

		
	
open FILE, "/home/ephas/Group/KOREAN/SF.csv" or die $!;
print "file opened\n";



my %raw; # hash of hash of hash: 1st key (person), 2nd key (model), 3rd key (protein\tsnp) and value (raw score)
my $i=0;
while (<FILE>) {
	chomp;
	unless ($_ =~ /^Exome/){
		my @tmp = split(/\t/,$_);
		$raw{$tmp[0]}{$tmp[3]}{"$tmp[1]\t$tmp[2]"}=$tmp[4];
		
		#last if(!(++$i%10000));
	}
}

# obs
print "obs\n";
my $llr = calc(\%raw,$pheno);
my %llr = %{$llr}; 
#print Dumper \%llr;

# randomized
print "randomized\n";
my %raw_shuffle;# hash of hash of hash: 1st key (person), 2nd key (model), 3rd key (protein\tsnp) and value (shuffled raw score)
%raw_shuffle=shuffle(\%raw);
my $llr_shuffle = calc(\%raw_shuffle,$pheno);
my %llr_shuffle = %{$llr_shuffle}; 
print Dumper \%llr_shuffle;

sub shuffle{
	my %raw = %{(shift)}; # hash of hash of hash: 1st key (person), 2nd key (model), 3rd key (protein\tsnp) and value (raw score)

	my %raw_shuffle;# hash of hash of hash: 1st key (person), 2nd key (model), 3rd key (protein\tsnp) and value (shuffled raw score)

	my %raw_index;# hash of hash of hash: 1st key (person), 2nd key (model), 3rd key (protein\tsnp) and value (index)
	my @index_score;
	my $i=0;
	foreach my $person (keys %raw){
		foreach my $model (keys %{$raw{$person}}){
			foreach my $snp (keys %{$raw{$person}{$model}}){
				$raw_index{$person}{$model}{$snp}=$i;
				push @index_score,$raw{$person}{$model}{$snp};
				$i++;
			}
		}
	}

	fisher_yates_shuffle( \@index_score );
	
	foreach my $person (keys %raw){
		foreach my $model (keys %{$raw{$person}}){
			foreach my $snp (keys %{$raw{$person}{$model}}){
				$raw_shuffle{$person}{$model}{$snp}=$index_score[$raw_index{$person}{$model}{$snp}];
			}
		}
	}
	
	return %raw_shuffle;
}

sub fisher_yates_shuffle {
    my $array = shift;
    my $i;
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
    }
    
}


sub calc{

	my %raw = %{(shift)}; # hash of hash of hash: 1st key (person), 2nd key (model), 3rd key (protein\tsnp) and value (raw score)
	my $pheno = shift;
	
	my %bscore_pos;
	my %bscore_neg;
	my %min_pos; 
	my %max_pos; 
	my %min_neg; 
	my %max_neg;
	my %PO_count_pos;
	my %PO_count_neg;
	my %count_pos;
	my %count_neg;
	my	$dbh = rackham::DBConnect('superfamily');
	
	my %llr;
	
	foreach my $person (keys %raw){	    	
			
		unless(exists($min_neg{$person})){
				$min_pos{$person} = 999999999; 
				$max_pos{$person} = 0; 
				$min_neg{$person}= 999999999; 
				$max_neg{$person} = 0;
				$count_pos{$person} = 0;
				$count_neg{$person} = 0;
		}
		
		foreach my $model (keys %{$raw{$person}}){
			my (  $sth );
			#print "Before query\n";
			$sth =   $dbh->prepare( "select distinct model.model,PO_mapping.po,all_score from PO_mapping,model where model.sf = PO_mapping.id and PO_mapping.inherited_from !='' and model.model=\"$model\" and PO_mapping.obo in ('$pheno');" );
		    $sth->execute;
			#print "After query\n";
			while (my ($sf,$PO,$all_score)= $sth->fetchrow_array ) {
				foreach my $snp (keys %{$raw{$person}{$model}}){
					my $c = $snp;
					if($raw{$person}{$model}{$snp}<1){
						my $t = log(1/$raw{$person}{$model}{$snp})/log(2);
						
						$bscore_pos{$person}{$PO}{$c} = $t;
						if ($t < $min_pos{$person}){
							$min_pos{$person} = $t;
						}
						if ($t > $max_pos{$person}){
							$max_pos{$person} = $t;
						}
						if(exists($PO_count_pos{$person}{$PO})){
							$PO_count_pos{$person}{$PO} = $PO_count_pos{$person}{$PO}+1
						}else{
							$PO_count_pos{$person}{$PO} = 1;
						}
						$count_pos{$person}++;
					}else{
						my $t = log($raw{$person}{$model}{$snp})/log(2);
						$bscore_neg{$person}{$PO}{$c} = $t;
						if ($t < $min_neg{$person}){
							$min_neg{$person} = $t;
						}
						if ($t > $max_neg{$person}){
							$max_neg{$person} = $t;
						}
						if(exists($PO_count_neg{$person}{$PO})){
							$PO_count_neg{$person}{$PO} = $PO_count_neg{$person}{$PO}+1
						}else{
							$PO_count_neg{$person}{$PO} = 1;
						}
						$count_neg{$person}++;
					}	
					
			    }				
			}

		}
		
		
		my %pos;
		foreach my $person (keys %bscore_pos){
			foreach my $PO (keys %{$bscore_pos{$person}}){
				foreach my $SNP (keys %{$bscore_pos{$person}{$PO}}){
					my $ascore = ($PO_count_pos{$person}{$PO}/$count_pos{$person})*(($bscore_pos{$person}{$PO}{$SNP}-$min_pos{$person})/($max_pos{$person}-$min_pos{$person}));
					if(exists($pos{$person}{$PO})){
					$pos{$person}{$PO} = $pos{$person}{$PO}+$ascore;
					}else{
						$pos{$person}{$PO} = $ascore;
					}
					
				}
			}
		}
		my %neg;
		foreach my $person (keys %bscore_neg){
			foreach my $PO (keys %{$bscore_neg{$person}}){
				foreach my $SNP (keys %{$bscore_neg{$person}{$PO}}){
					my $ascore = ($PO_count_neg{$person}{$PO}/$count_neg{$person})*(($bscore_neg{$person}{$PO}{$SNP}-$min_neg{$person})/($max_neg{$person}-$min_neg{$person}));
					if(exists($neg{$person}{$PO})){
					$neg{$person}{$PO} = $neg{$person}{$PO}+$ascore;
					}else{
						$neg{$person}{$PO} = $ascore;
					}
				}
			}
		}
		my %prob;
		foreach my $person (%pos){
			foreach my $PO (keys %{$pos{$person}}){
					my $pos_score = $pos{$person}{$PO};
					my $neg_score;
					if(exists($neg{$person}{$PO})){
						$neg_score = $neg{$person}{$PO};
						$prob{$person}{$PO} = $pos_score/($pos_score+$neg_score);
					}else{
						$prob{$person}{$PO} = 1;
					}
					my $llr;
					if($prob{$person}{$PO} == 1){
					$llr = 100000;	
					$llr{$person}{$PO} = $llr;
					}else{
					$llr = $prob{$person}{$PO}/(1-$prob{$person}{$PO});
					$llr{$person}{$PO} = $llr;
					}
			}
		}
		
		
	}
	return \%llr;
}
