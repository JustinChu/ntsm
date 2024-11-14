#!/usr/bin/perl
use warnings;
use strict;
use diagnostics;
use IO::File;

my $fh     = new IO::File( $ARGV[0], "r" );
my $prefix = $ARGV[1];
my $line   = $fh->getline();

#id->type(AT,CG)->count
my %idToUniqCount;
my %idToStr;

my $maxCount = 1;

while ($line) {
	chomp($line);
	if ($line) {

#rs4970383|0|AT	0	chr1	903160	20	19M	*	0	0	CCCATGACCCCATTTCACC	*	XT:A:U	NM:i:0	X0:i:1	X1:i:2	XM:i:0	XO:i:0	XG:i:0	MD:Z:19	XA:Z:chr20,+44466325,19M,1;chr15,+45228912,19M,1;
		my @tmpAr = split( /\t/, $line );
		if ( $tmpAr[0] =~ /([^\|]+)\|\d+\|(AT|CG)/ ) {
			my $id   = $1;
			my $type = $2;
			my $seq  = $tmpAr[9];
			if ( !exists( $idToUniqCount{$id}{$type} ) ) {
				$idToUniqCount{$id}{$type} = 0;
			}
			if ( $line =~ /X0:i:(\d+)/ ) {
				my $count = $1;
				if ( $line =~ /.*X1:i:(\d+)/ ) {
					$count += $1;
				}
				if ( $count > 1 ) {
					if ( $maxCount < ++$idToUniqCount{$id}{$type} ) {
						$maxCount = $idToUniqCount{$id}{$type};
					}
				}
				else {
					if ( exists( $idToStr{$id}{$type} ) ) {
						$idToStr{$id}{$type} .= 'N' . $seq;
					}
					else {
						$idToStr{$id}{$type} = $seq;
					}
				}
			}
			else {
				if ( exists( $idToStr{$id}{$type} ) ) {
					$idToStr{$id}{$type} .= 'N' . $seq;
				}
				else {
					$idToStr{$id}{$type} = $seq;
				}
			}
		}
		else {
			print STDERR "unable to parse: " . $line . "\n";
		}
	}
	$line = $fh->getline();
}
$fh->close();

my @fhs;

for ( my $i = 0 ; $i < $maxCount ; ++$i ) {
	push( @fhs, new IO::File( $prefix . "_n" . $i . ".fa", 'w' ) );
}

foreach my $id ( sort keys(%idToUniqCount) ) {
	for ( my $i = 0 ; $i < $maxCount ; ++$i ) {
		if (   $idToUniqCount{$id}{"AT"} <= $i
			&& $idToUniqCount{$id}{"CG"} <= $i )
		{
			if (   exists( $idToStr{$id}{"AT"} )
				&& exists( $idToStr{$id}{"CG"} ) )
			{
				$fhs[$i]->write( ">" . $id . " ref\n" . $idToStr{$id}{"AT"} . "\n" );
				$fhs[$i]->write( ">" . $id . " var\n" . $idToStr{$id}{"CG"} . "\n" );
			}
			else {
				print STDERR $id . " " . $i . "\n";
			}
		}
	}
}
