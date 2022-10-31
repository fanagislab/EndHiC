#!/usr/bin/perl
use strict;

# parse input options
if (not $ARGV[0] or $ARGV[0] eq "-h") {
	die "
	Description: read gfa.cluster produced by endhic_pipeline.pl, select the final clusters as chromosomes, and output the ordered and oriented contigs in AGP format.
	
	Author: Sen Wang, wangsen1993@163.com, 2021/7/26.

	Usage: perl cluster2agp.pl gfa.cluster contigs.len > gfa.cluster.agp
	\n";
}

# read gfa.cluster
my %chr;
my ($num_clu, $len_clu) = (0, 0);
open IN, "<$ARGV[0]" or die "Cannot open $ARGV[0]!\n";
while (<IN>) {
	next if /^#/;
	chomp;
	my @f = split(/[\t;]/, $_);
	$num_clu += $f[1];
	$len_clu += $f[2];
	@{$chr{$f[0]}} = @f[4 .. @f - 1];
}
close IN;

# read contigs.len
my %lens;
my ($num_tot, $len_tot) = (0, 0);
open IN, "<$ARGV[1]" or die "Cannot open $ARGV[1]!\n";
while (<IN>) {
	chomp;
	my ($c, $l) = split(/\t/, $_);
	$lens{$c} = $l;
	$num_tot += 1;
	$len_tot += $l;
}
close IN;

# output AGP
my $number = keys %chr;
foreach my $clu (sort keys %chr) {
	my $start = 1;
	my $end = 0;
	my $rank = 0;
	my $ctg = "";
	my $strand = "";
	my @tem = @{$chr{$clu}};
	if (@tem == 1) {
		$rank += 1;
		$ctg = substr($tem[0], 0, length($tem[0]) - 1);
		die "Cannot get the length of $ctg! check $ARGV[1]!\n" if not $lens{$ctg};
		$end = $start + $lens{$ctg} - 1;
		$strand = substr($tem[0], -1, 1);
		$clu =~ s/[A-Za-z]+_[A-Z]?/Chr_/;
		print "$clu\t$start\t$end\t$rank\tW\t$ctg\t1\t$lens{$ctg}\t$strand\n";
		delete $lens{$ctg};
	} else {
		foreach my $c (@tem[0 .. @tem - 2]) {
			$ctg = substr($c, 0, length($c) - 1);
			$strand = substr($c, -1, 1);
			$rank += 1;
			die "Cannot get the length of $ctg! check $ARGV[1]!\n" if not $lens{$ctg};
			$end = $start + $lens{$ctg} - 1;
			$clu =~ s/[A-Za-z]+_[A-Z]?/Chr_/;
			print "$clu\t$start\t$end\t$rank\tW\t$ctg\t1\t$lens{$ctg}\t$strand\n";
			delete $lens{$ctg};
			$start = $end + 1;
			$end = $start + 1000 - 1;
			$rank += 1;
			print "$clu\t$start\t$end\t$rank\tU\t1000\tcontig\tyes\tmap\n";
			$start = $end + 1;
		}
		$ctg = substr($tem[-1], 0, length($tem[-1]) - 1);
		die "Cannot get the length of $ctg! check $ARGV[1]!\n" if not $lens{$ctg};
		$strand = substr($tem[-1], -1, 1);
		$rank += 1;
		$end = $start + $lens{$ctg} - 1;
		print "$clu\t$start\t$end\t$rank\tW\t$ctg\t1\t$lens{$ctg}\t$strand\n";
		delete $lens{$ctg};
	}
}
foreach my $c (sort keys %lens) {
	print "$c\t1\t$lens{$c}\t1\tW\t$c\t1\t$lens{$c}\t\+\n";
}

# output length stat
my $rate = sprintf("%.2f", $len_clu / $len_tot * 100);
print STDERR "$num_clu of $num_tot contigs were anchored into $number chromosomes.\n";
print STDERR "Total chromosome length $len_clu accounted for $rate\% of total contig length $len_tot.\n";
