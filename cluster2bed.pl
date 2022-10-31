#!/usr/bin/perl
use strict;
if (not $ARGV[0] or $ARGV[0] eq "-h") {
	die "
	Description: read gfa.cluster and matrix_abs.bed file, reorder bins according to the ordered and oriented contigs in each cluster, and output new cluster.bed and cluster.len.

	Author: Sen Wang, wangsen1993@163.com, 2021/9/2.

	Usage: perl cluster2bed.pl formal_100000_abs.bed z.EndHiC.results.summary.and.analysis.A.cluster >z.EndHiC.results.summary.and.analysis.A.cluster_abs.bed 2>z.EndHiC.results.summary.and.analysis.A.cluster.id.len
	\n";
}
# read bed and map bins to contigs
my %ctg2bin;
open IN, "<$ARGV[0]" or die "Cannot open $ARGV[0]!\n";
while (<IN>) {
	chomp;
	my @f = split(/\t/, $_);
	my $c = shift @f;
	push @{$ctg2bin{$c}}, join("\t", @f);
}
close IN;
# read cluster and map contigs to clusters
my %clu2ctg;
open IN, "<$ARGV[1]" or die "Cannot open $ARGV[1]!\n";
while (<IN>) {
	next if /^#/;
	chomp;
	my @f = split(/\t/, $_);
	print STDERR "$f[0]\t$f[2]\n";
	$clu2ctg{$f[0]} = $f[4];
}
close IN;
# reorder the bins according to the mapped contigs in each cluster
foreach my $clu (sort keys %clu2ctg) {
	my @ctgs = split(/;/, $clu2ctg{$clu});
	my $start = 0;
	my $end = 0;
	#print STDOUT "@ctgs\n";
	foreach my $c (@ctgs) {
		my $ctg = substr($c, 0, length($c) - 1);
		die "Cannot get the bins of $ctg! check $ARGV[0]!\n" if not $ctg2bin{$ctg};
		my $strand = substr($c, -1, 1);
		if ($strand eq "+") {
			foreach my $line (@{$ctg2bin{$ctg}}) {
				my ($s, $e, $i) = split(/\t/, $line);
				$end = $start + ($e - $s);
				print STDOUT "$clu\t$start\t$end\t$i\n";
				$start = $end;
			}
		} elsif ($strand eq "-") {
			foreach my $line (reverse @{$ctg2bin{$ctg}}) {
				my ($s, $e, $i) = split(/\t/, $line);
				$end = $start + ($e - $s);
				print STDOUT "$clu\t$start\t$end\t$i\n";
				$start = $end;
			}
		}
	}
}
