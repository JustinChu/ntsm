#!/usr/bin/perl
use warnings;
use strict;
use diagnostics;
use IO::File;

#rsID -> chr pos
my %chrPosToRSID;

#vcf file with IDs to filter main list with
my $fh   = new IO::File( $ARGV[0], "r" );
my $line = $fh->getline();

while ($line) {
	if ($line) {
		chomp($line);

		#chr1	918870	rs7537756	A	G	.	.
		my @tempArr = split( /\s/, $line );
		$chrPosToRSID{ $tempArr[0] . "\t" . $tempArr[1] } = $tempArr[2];
	}
	$line = $fh->getline();
}
$fh->close();

#genomeIDs file
my $fh2 = new IO::File( $ARGV[1], "r" );
$line = $fh2->getline();

while ($line) {
	print "rsID\t" . $line . "\n";
	$line = $fh2->getline();
}
$fh2->close();

#main list VCF
my $fh3 = new IO::File( $ARGV[1], "r" );
$line = $fh3->getline();

my $count = 0;

while ($line) {
	if ($line) {
		chomp($line);
		my @tempArr = split( /\s/, $line );
		my $id      = $tempArr[0] . "\t" . $tempArr[1];
		if ( exists( $chrPosToRSID{$id} ) ) {
			print $chrPosToRSID{$id};
			for ( my $i = 9 ; $i < scalar(@tempArr) ; ++$i ) {
				print "\t". $tempArr[$1];
			}
			print "\n";
			$count++;
		}
	}
	$line = $fh3->getline();
}
$fh3->close();

print STDERR $count . " vs " . scalar(%chrPosToRSID) . "\n";

