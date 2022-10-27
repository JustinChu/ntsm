#!/usr/bin/perl
use warnings;
use strict;
use diagnostics;
use IO::File;

#rsID -> chr pos
my %chrPosToRSID;
my %chrPosToVar;

#vcf file with IDs to filter main list with
my $fh   = new IO::File( $ARGV[0], "r" );
my $line = $fh->getline();

while ($line) {
	if ($line) {
		chomp($line);

		#chr1	918870	rs7537756	A	G	.	.
		my @tempArr = split( /\s/, $line );
		$chrPosToRSID{ $tempArr[0] . "\t" . $tempArr[1] } = $tempArr[2];
		$chrPosToVar{ $tempArr[0] . "\t" . $tempArr[1] }  = $tempArr[4];
	}
	$line = $fh->getline();
}
$fh->close();

#genomeIDs file
my $fh2 = new IO::File( $ARGV[1], "r" );
$line = $fh2->getline();

while ($line) {
	chomp($line);
	print "rsID\t" . $line . "\n";
	$line = $fh2->getline();
}
$fh2->close();

#main list VCF
my $fh3 = new IO::File( $ARGV[2], "r" );
$line = $fh3->getline();
my $count = 0;

while ($line) {
	if ($line) {
		chomp($line);

		#chr1	911428	1:911428:C:T	C	T
		my @tempArr = split( /\t/, $line );
		my $id      = $tempArr[0] . "\t" . $tempArr[1];
		my $refBase = $tempArr[3];
		my $varBase = $tempArr[4];
		if ( exists( $chrPosToRSID{$id} ) ) {
			if ( $chrPosToVar{$id} eq $varBase ) {
				print $chrPosToRSID{$id};
				if ( $refBase eq "A" || $refBase eq "T" ) {
					for ( my $i = 9 ; $i < scalar(@tempArr) ; ++$i ) {
						if ( $tempArr[$i] eq "1|1" ) {
							print "\t2";
						}
						elsif ( $tempArr[$i] eq "1|0" || $tempArr[$i] eq "0|1" )
						{
							print "\t1";
						}
						elsif ( $tempArr[$i] eq "0|0" ) {
							print "\t0";
						}
						else {
							print "\tNA";
						}
					}
				}
				else {
					for ( my $i = 9 ; $i < scalar(@tempArr) ; ++$i ) {
						if ( $tempArr[$i] eq "0|0" ) {
							print "\t2";
						}
						elsif ( $tempArr[$i] eq "1|0" || $tempArr[$i] eq "0|1" )
						{
							print "\t1";
						}
						elsif ( $tempArr[$i] eq "1|1" ) {
							print "\t0";
						}
						else {
							print "\tNA";
						}
					}
				}
				print "\n";
				$count++;
			}
		}
	}
	$line = $fh3->getline();
}
$fh3->close();

print STDERR $count . " vs " . scalar(%chrPosToRSID) . "\n";

