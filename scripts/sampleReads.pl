#!/usr/bin/perl

#given reference genomeID and file -> extract from master bed file

use warnings;
use strict;
use diagnostics;
use IO::File;

my $fh         = new IO::File( $ARGV[0], "r" );
my $maxCount   = 6000000000;
my $curHeader  = $fh->getline();
my $curSeq     = $fh->getline();
my $spacer     = $fh->getline();
my $curQual    = $fh->getline();
my $totalBases = 0;
while ( $maxCount > $totalBases ) {
	print $curHeader . $curSeq . $spacer . $curQual;
	$totalBases += length($curSeq)-1;
	$curHeader = $fh->getline();
	$curSeq    = $fh->getline();
	$spacer    = $fh->getline();
	$curQual   = $fh->getline();
}
$fh->close();
