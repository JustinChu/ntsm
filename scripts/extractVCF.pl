#!/usr/bin/perl

#given reference genomeID and file -> extract from master bed file

use warnings;
use strict;
use diagnostics;
use IO::File;

#CHM13#0#chr6	31786100	31911587	*	0	+
my $fh   = new IO::File( $ARGV[0], "r" );
my $line = $fh->getline();

while ($line) {
	chomp($line);
	my @tempArr = split( /\t/, $line );

#1865	chr1	167880175	167880176	rs203849	0	-	A	A	C/T	genomic	single	by-cluster,by-frequency,by-submitter,by-2hit-2allele,by-1000genomes	0.499908	0.006788	coding-synon	exact	1	InconsistentAlleles	25	1000GENOMES,BCM-HGSC-SUB,BCMHGSC_JDW,BUSHMAN,CLINSEQ_SNP,COMPLETE_GENOMICS,CORNELL,ENSEMBL,EVA-GONL,EXOME_CHIP,GMI,HUMANGENOME_JCVI,ILLUMINA,ILLUMINA-UK,JMKIDD_LAB,KRIBB_YJKIM,KWOK,NHLBI-ESP,PERLEGEN,PJP,SC_JCM,SEATTLESEQ,SSMP,TISHKOFF,TSC-CSHL,	0				maf-5-some-pop,genotype-conflict

	if ( $tempArr[6] eq "-" ) {
		my $secondaryAllele = revComp( $tempArr[9] );
		$secondaryAllele =~ s/\Q$tempArr[7]//g;
		$secondaryAllele =~ s/\///g;
		print $tempArr[1] . "\t"
		  . $tempArr[3] . "\t"
		  . $tempArr[4] . "\t"
		  . $tempArr[7] . "\t"
		  . $secondaryAllele . "\n";
	}
	else {
		my $secondaryAllele = $tempArr[9];
		$secondaryAllele =~ s/\Q$tempArr[7]//g;
		$secondaryAllele =~ s/\///g;
		print $tempArr[1] . "\t"
		  . $tempArr[3] . "\t"
		  . $tempArr[4] . "\t"
		  . $tempArr[7] . "\t"
		  . $secondaryAllele . "\n";
	}

	$line = $fh->getline();
}
$fh->close();

sub revComp {
	my $revcomp = reverse shift;
	$revcomp =~ tr/ATGCatgc/TACGtacg/;
	return $revcomp;
}
