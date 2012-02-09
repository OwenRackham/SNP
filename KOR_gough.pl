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
my %POs;
open DIST,">dist.txt";
while (<FILE>) {
	    
	    
		unless ($_ =~ /^Exome/){
		my @line = split(/\t/,$_);
		my $t = -log($line[4])/log(2);
		print DIST "$t\n";
		my (  $sth );
		#print "Before query\n";
		$sth =   $dbh->prepare( "select distinct model.sf,PO_mapping.po,all_score from PO_mapping,model where model.sf = PO_mapping.id and PO_mapping.inherited_from !='' and model.model = $line[3] and PO_mapping.obo in ('$pheno');" );
        #$sth->execute;
		#print "After query\n";
		my $i=0;
        #while (my ($sf,$PO,$all_score)= $sth->fetchrow_array ) {
        	
        #}		
		print "$line[3]:\t$i\n";
	}
}

