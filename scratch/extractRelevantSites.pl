#!/usr/bin/perl
use warnings;
use strict;
use diagnostics;
use IO::File;

#rsID -> chr pos
my %chrPosToRSID;

#vcf file with IDs to filter main list with
my $fh = new IO::File( $ARGV[0], "r" );
my $line = $fh->getline();

while ($line) {
	if ($line) {
		chomp($line);

		#chr1	918870	rs7537756	A	G	.	.
		my @tempArr = split( /\s/, $line );
		$chrPosToRSID{ $tempArr[0] . "\t"
			  . $tempArr[1] . "\t"
			  . $tempArr[3] . "\t"
			  . $tempArr[4] } = $tempArr[2];
	}
	$line = $fh->getline();
}
$fh->close();

#main list VCF
my $fh2 = new IO::File( $ARGV[1], "r" );
$line = $fh2->getline();
my $headerStr = "";



while ($line) {
	if ($line) {
		chomp($line);
		if ( $line !~ /^#/ ) {
			#chr1	911428	1:911428:C:T	C	T
			my @tempArr = split( /\t/, $line );
			my $id =
			    $tempArr[0] . "\t"
			  . $tempArr[1] . "\t"
			  . $tempArr[3] . "\t"
			  . $tempArr[4];
			if ( exists( $chrPosToRSID{$id} ) ) {
				$tempArr[2] = $chrPosToRSID{$id};
				print $tempArr[0];
				for ( my $i = 1 ; $i < scalar(@tempArr) ; ++$i ) {
					print "\t" . $tempArr[$i];
				}
				print "\n";
			}
		}
		elsif($headerStr eq ""){
			if($line =~ /^#CHROM/){
				$headerStr = $line;
				print $line ."\n";
			}
		}
	}
	$line = $fh2->getline();
}
$fh2->close();

