#!/usr/bin/perl
use warnings;
use strict;
use diagnostics;
use IO::File;

my $fh       = new IO::File( $ARGV[0], "r" );
my $prefix   = $ARGV[1];
my $maxCount = $ARGV[2] - $ARGV[3] + 1;
my $line     = $fh->getline();
chomp($line);

#id->type(AT,CG)->count
my %idToUniqCount;
my %idToStr;

#my %currPos = (
#	"AT" => 0,
#	"CG" => 0
#);

while ($line) {
	chomp($line);

#rs4970383|0|AT	0	chr1	903160	20	19M	*	0	0	CCCATGACCCCATTTCACC	*	XT:A:U	NM:i:0	X0:i:1	X1:i:2	XM:i:0	XO:i:0	XG:i:0	MD:Z:19	XA:Z:chr20,+44466325,19M,1;chr15,+45228912,19M,1;
	my @tmpAr = split( /\t/, $line );
	if ( $tmpAr[0] =~ /([^\|]+)\|(\d+)\|(AT|CG)/ ) {
		my $id   = $1;
		my $pos  = $2;
		my $type = $3;
		my $seq  = $tmpAr[9];
		if ( !exists( $idToUniqCount{$id}{$type} ) ) {
			$idToUniqCount{$id}{$type} = $maxCount;
		}
		if ( $line =~ /X0:i:(\d+)/ ) {
			my $count = $1;
			if ( $line =~ /.*X1:i:(\d+)/ ) {
				$count += $1;
			}
			if ( $count > 1 ) {
				if ( exists( $idToStr{$id}{$type} ) ) {
#					if ( $currPos{$type} - $pos == 1 ) {
#						$idToStr{$id}{$type} .= substr( $seq, -1 );
#					}
#					else {
						$idToStr{$id}{$type} .= 'N' . $seq;
#					}
				}
				else {
					$idToStr{$id}{$type} = $seq;
				}
#				$currPos{$type} = $pos;
				--$idToUniqCount{$id}{$type};
			}
		}
		else {
			if ( exists( $idToStr{$id}{$type} ) ) {
#				if ( $currPos{$type} - $pos == 1 ) {
#					$idToStr{$id}{$type} .= substr( $seq, -1 );
#				}
#				else {
					$idToStr{$id}{$type} .= 'N' . $seq;
#				}
			}
			else {
				$idToStr{$id}{$type} = $seq;
			}
#			$currPos{$type} = $pos;
			--$idToUniqCount{$id}{$type};
		}
	}
	else {
		print STDERR "unable to parse: " . $line . "\n";
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
		if (   exists( $idToUniqCount{$id}{"AT"} )
			&& exists( $idToUniqCount{$id}{"CG"} ) )
		{
			if (   $idToUniqCount{$id}{"AT"} <= $i
				&& $idToUniqCount{$id}{"CG"} <= $i )
			{
				if (   exists( $idToStr{$id}{"AT"} )
					&& exists( $idToStr{$id}{"CG"} ) )
				{
					$fhs[$i]->write(
						">" . $id . " ref\n" . $idToStr{$id}{"AT"} . "\n" );
					$fhs[$i]->write(
						">" . $id . " var\n" . $idToStr{$id}{"CG"} . "\n" );
				}
				else {
					print STDERR "Possible file truncation. Missing: "
					  . $id . " "
					  . $i . " "
					  . $idToUniqCount{$id}{"AT"} . " "
					  . $idToUniqCount{$id}{"CG"} . " "
					  . int( exists( $idToStr{$id}{"AT"} ) ) . " "
					  . int( exists( $idToStr{$id}{"CG"} ) ) . "\n";
				}
			}
		}
	}
}
