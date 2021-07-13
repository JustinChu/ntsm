#!/usr/bin/perl

#given reference genomeID and file -> extract from master bed file

use warnings;
use strict;
use diagnostics;
use IO::File;

my $fh   = new IO::File( $ARGV[0], "r" );
my $line = $fh->getline();

my $threshold = 3;

while ($line) {
	chomp($line);
	my @tempArr = split( /\s/, $line );
	
	if($tempArr[0] > $threshold && $tempArr[2] > $threshold){
		print "3";
	}
	elsif($tempArr[0] > $threshold){
		print "1";
	}
	elsif($tempArr[1] > $threshold){
		print "2";
	}
	else{
		print "0";
	}
	print "\n";

	$line = $fh->getline();
}
$fh->close();