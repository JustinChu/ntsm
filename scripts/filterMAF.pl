#!/usr/bin/perl
#Filters out MAF < 0.05, keeps only the 2 highest alleles per site (Reorders only the 4th and 5th columns)
#assumes dbGaP_PopFreq is present on each entry, otherwise it ignores the line

use warnings;
use strict;
use diagnostics;
use IO::File;

my $fh        = new IO::File( $ARGV[0], "r" );
my $line      = $fh->getline();
my $threshold = 0.05;

my $sum   = 0;
my $count = 0;

while ($line) {
	chomp($line);
	if ($line) {    #dbGaP_PopFreq:0.07241,0.9276,0;COMMON
		if ( $line =~ /dbGaP_PopFreq:([^;]+)/ ) {
			my $alleleStr = $1;
			my @alleleArr = split( /,/,  $alleleStr );
			my @lineArr   = split( /\t/, $line );
			my @ntArr     = ();
			push( @ntArr, $lineArr[3] );
			push( @ntArr, split( /,/, $lineArr[4] ) );

			my $highestAllele       = 0;
			my $secondHighestAllele = 0;
			my $highestNT           = $ntArr[0];
			my $secondNT            = $ntArr[1];

			for ( my $i = 0 ; $i < scalar(@alleleArr) ; ++$i ) {
				my $allele = $alleleArr[$i];
				unless ( $allele eq '.' ) {
					if ( $allele > $highestAllele ) {
						$secondHighestAllele = $highestAllele;
						$highestAllele       = $allele;
						$secondNT            = $highestNT;
						$highestNT           = $ntArr[$i];
					}
					elsif ( $allele > $secondHighestAllele ) {
						$secondHighestAllele = $highestAllele;
						$secondNT            = $ntArr[$i];
					}
				}
			}
			$lineArr[3] = $highestNT;
			$lineArr[4] = $secondNT;

			if ( $secondHighestAllele > $threshold ) {
				local $" = "\t";
				print "@lineArr\n";
			}
		}
		else {
			print STDERR $line . "\n";
		}
	}
	$count++;
	$line = $fh->getline();
}
$fh->close();
