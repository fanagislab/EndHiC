#!/usr/bin/perl
use strict;
use Getopt::Long;

# parse input options
my ($help, $unclu);
GetOptions (
	'un|unclustered' => \$unclu,
	'h|help!' => \$help
);
if (not $ARGV[0] or $help) {
	die "
	Description: read GFA clusters at contig and cluster level and output the new connected clusters by their belonging contigs.

	Author: Sen Wang, wangsen1993@163.com, 2021/8/16

	Usage: perl cluster_merge.pl [--unclustered] contigs_gfa.cluster cluster_gfa.cluster > merged_gfa.cluster
	\n";
}

# read GFA cluster at contig level
my (%cluster, %num, %len, %robust);
open IN, "<$ARGV[0]" or die "Cannot open $ARGV[0]!\n";
while (<IN>) {
	chomp;
	next if /^#/;
	my @f = split(/\t/, $_);
	my @l = split(/;/, $f[4]);
	@{$cluster{$f[0]}} = @l;
	$num{$f[0]} = $f[1];
	$len{$f[0]} = $f[2];
	$robust{$f[0]} = $f[3];
}
close IN;

# read GFA cluster at cluster level and output the connected clusters
open IN, "<$ARGV[1]" or die "Cannot open $ARGV[1]!\n";
while (<IN>) {
	chomp;
	if (/^#/) {
		print "$_\n";
		next;
	}
	my @f = split(/\t/, $_);
	my @l = split(/;/, $f[4]);
	if (@l > 1) {
		my @tem;
		foreach my $c (@l) {
			if ($c =~ /\+$/) {
				$c =~ s/[\+\-]//;
				push @tem, @{$cluster{$c}};
				delete $cluster{$c};
			} elsif ($c =~ /\-$/) {
				my @rev;
				$c =~ s/[\+\-]//;
				foreach my $ctg (reverse @{$cluster{$c}}) {
					if ($ctg =~ /\+$/) {
						$ctg =~ s/\+/\-/;
						push @rev, $ctg;
					} elsif ($ctg =~ /\-$/) {
						$ctg =~ s/\-/\+/;
						push @rev, $ctg;
					}
				}					
				push @tem, @rev;
				delete $cluster{$c};
			}
		}
		$f[4] = join(";", @tem);
		$f[1] = @tem;
		print join("\t", @f) . "\n";
	} elsif (@l == 1) {
		$f[4] =~ s/[\+\-]//;
		my $k = $f[4];
		my @tem = @{$cluster{$k}};
		$f[1] = @tem;
		$f[4] = join(";", @tem);
		print join("\t", @f) . "\n";
		delete $cluster{$k};
	}
}
close IN;

# output the unconnected clusters
my @n = keys %cluster;
if ($unclu and @n > 0) {
	foreach my $k (sort {$len{$b} <=> $len{$a}} keys %cluster) {
		print "$k\t$num{$k}\t$len{$k}\t$robust{$k}\t" . join(";", @{$cluster{$k}}) . "\n";
	}
}
