#!/usr/bin/perl

use warnings;
use strict;
use diagnostics;
use IO::File;

my @dataStr;
my @countA;
my @countB;

while (<>) {
	my $fh    = new IO::File( $_, "r" );
	my $line  = $fh->getline();
	my $index = 0;
	while ($line) {
		chomp($line);
		if ($line) {

			#rs431480	1	24
			my @tmpArr = split( /\t/, $line );
			if ( scalar(@dataStr) == $index ) {
				push( @dataStr, $tmpArr[0] );
				push( @countA,  $tmpArr[1] );
				push( @countB,  $tmpArr[2] );
			}
			else {
				$dataStr[$index] += $tmpArr[0];
				$countA[$index]  += $tmpArr[1];
				$countB[$index]  += $tmpArr[2];
			}
		}
		$index++;
		$line = $fh->getline();
	}
	$fh->close();
}

for ( my $i = 0 ; $i < scalar(@dataStr) ; ++$i ) {
	print $dataStr[$i] . "\t"
	  . $countA[$i] . "\t"
	  . $countB[$i] + "\n";
}
