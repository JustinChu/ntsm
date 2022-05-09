#!/usr/bin/perl

#given reference genomeID and file -> extract from master bed file

use warnings;
use strict;
use diagnostics;
use IO::File;

my $fh        = new IO::File( $ARGV[0], "r" );
my $line      = $fh->getline();
my $threshold = 3;

my $sum = 0;
my $count = 0;

while ($line) {
	chomp($line);
	if ($line) {
		my @tempArr = split( /\s/, $line );
#		print $tempArr[1];
		$sum += $tempArr[1] + $tempArr[2];

		if ( $tempArr[1] > $threshold && $tempArr[2] > $threshold ) {
			print "3";	
		}
		elsif ( $tempArr[1] > $threshold ) {
			print "1";
		}
		elsif ( $tempArr[2] > $threshold ) {
			print "2";
		}
		else {
			print "0";
		}
	}
	$count++;
	$line = $fh->getline();
}
print "\n";

my $meanCov = ($sum/($count*2));
print $meanCov. "\n";
$fh->close();
