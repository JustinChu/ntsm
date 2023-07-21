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
my $mafSum = 0;
my $mafCount = 0;

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
			my $highestNT           = $lineArr[3];
			my $secondNT            = $ntArr[0];

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
						$secondHighestAllele = $allele;
						$secondNT            = $ntArr[$i];
					}
				}
			}
			if($lineArr[3] eq $highestNT){
				$lineArr[4] = $secondNT;
			}
			else{
#				unless($lineArr[3] eq $secondNT){
#					print STDERR "Reference Allele is not in high abundance:" . $line . "\n";
#				}
				$lineArr[4] = $highestNT;
			}

			if ( $secondHighestAllele > $threshold ) {
#				if($secondHighestAllele > 0.5){
#					print STDERR "Allele Frequency error: " . $line . "\n";
#				}
				++$mafCount;
				$mafSum += $secondHighestAllele;
				local $" = "\t";
				print "@lineArr\n";
			}
			else{
				$count++;
			}
		}
		else {
#			print STDERR 'Unparsed Line: ' . $line . "\n";
		}
	}
	$line = $fh->getline();
}
$fh->close();
print STDERR 'Total Removed: ' . $count . "\n";
print STDERR 'Avg MAF: ' . ($mafSum/$mafCount) . "\n";
