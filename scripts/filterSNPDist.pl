#!/usr/bin/perl
use warnings;
use strict;
use diagnostics;
use IO::File;

my %chrPosToRSID;

#master list of common SNPs
my $fh       = new IO::File( $ARGV[0], "r" );
my $line     = $fh->getline();
my $distance = 31;

while ($line) {
	chomp($line);
	my @lineArr = split( /\t/, $line );
	$chrPosToRSID{ $lineArr[0] . "\t" . $lineArr[1] } = $lineArr[2];
	$line = $fh->getline();
}
$fh->close();

#main list VCF
my $fh2 = new IO::File( $ARGV[0], "r" );
$line = $fh2->getline();
my $count = 0;

#for each SNP confirm position
while ($line) {
	chomp($line);
	my @lineArr = split( /\t/, $line );
	my $start   = $lineArr[1] - $distance + 1;
	my $end     = $lineArr[1] + $distance;
	my $isGood  = 1;
	for ( my $pos = $start ; $pos < $end ; ++$pos ) {
		if ( $pos != $lineArr[1]
			&& exists( $chrPosToRSID{ $lineArr[0] . "\t" . $pos } ) )
		{
			$isGood = 0;
		}
	}
	if($isGood == 1){
		print $line;
	}
	$line = $fh2->getline();
}
$fh2->close();

