#!/usr/bin/perl
use strict;

# parse input options
if (not $ARGV[0] or $ARGV[0] eq "-h") {
	die "
	Description: read AGP produced by cluster2agp.pl and output the chromosome-level scaffold sequences in multi-line FASTA format.
	
	Author: Sen Wang, wangsen1993@163.com, 2021/7/26.

	Usage: perl agp2fasta.pl gfa.cluster.agp contigs.fasta > gfa.cluster.agp.fasta
	\n";
}

# read cluster.agp
my (%scaffold, %contig);
open IN, "<$ARGV[0]" or die "Cannot open $ARGV[0]!\n";
while (<IN>) {
	chomp;
	my @f = split(/\t/, $_);
	if ($f[0] ne $f[5]) {
		push @{$scaffold{$f[0]}}, "$f[5]$f[8]";
	} else {
		$contig{$f[0]} = $f[5];
	}
}
close IN;

# read contigs.fasta
my %seqs;
my $header;
open IN, "<$ARGV[1]" or die "Cannot open $ARGV[1]!\n";
while (<IN>) {
	chomp;
	if (/^>(\w+)/) {
		$header = $1;
	} else {
		$seqs{$header} .= $_;
	}
}
close IN;

# output chromosome-level scaffold sequences
foreach my $s (sort keys %scaffold) {
	print ">$s\n";
	my $seq = "";
	foreach my $c (@{$scaffold{$s}}) {
		my $ctg = substr($c, 0, length($c) - 1);
		my $strand = substr($c, -1, 1);
		if ($strand eq '+') {
			die "Cannot get the sequence of $ctg! check $ARGV[1]!\n" if not $seqs{$ctg};
			$seq .= $seqs{$ctg};
		} elsif ($strand eq '-') {
			die "Cannot get the sequence of $ctg! check $ARGV[1]!\n" if not $seqs{$ctg};
			my $tem = $seqs{$ctg};
			$tem = reverse($tem);
			$tem =~ tr/ATCG/TAGC/;
			$seq .= $tem;
		} else {
			$ctg =~ /(\d+)/;
			$seq .= "N" x $1;
		}
	}
	for (my $i = 0; $i < length($seq); $i += 60) {
		my $sub = substr($seq, $i, 60);
		print "$sub\n";
	}
}
foreach my $c (sort keys %contig) {
	print ">$c\n";
	my $seq = $seqs{$c};
	die "Cannot get the sequence of $c! check $ARGV[1]!\n" if not $seqs{$c};
	for (my $i = 0; $i < length($seq); $i += 60) {
		my $sub = substr($seq, $i, 60);
		print "$sub\n";
	}
}
