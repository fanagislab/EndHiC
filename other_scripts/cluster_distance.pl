#!/usr/bin/perl

=head1 Name

cluster_distance.pl  -- calculate the distance between clusters

=head1 Description

Using Hi-C links from all the included contigs of compared clusters

=head1 Version

  Author: Fan Wei, fanwei@caas.cn
  Version: 1.0,  Date: 2021/7/23
  Note:

=head1 Usage
  
  cluster_distance.pl  <endhic_result_file> <halfContig_contact_file>

  --help            output help information to screen  

=head1 Example

    perl cluster_distance.pl z.EndHiC.results.summary.and.analysis.B.cluster humanHiC_100000.matrix.100000.halfContig.CtgContact > z.EndHiC.results.summary.and.analysis.B.cluster.distance

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;


##get options from command line into variables and set default values
my ($Verbose,$Help);
GetOptions(
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if (@ARGV == 0 || $Help);


my $endhic_result_file = shift;
my $halfContig_contact_file = shift;

my %CtgCluster;
my %ClusterContact;
my %ClusterLen;

##Cluster_01      4       242687818       90      ptg000014l-;ptg000033l+;ptg000050l-;ptg000021l-
open IN,  $endhic_result_file || die "fail  $endhic_result_file";
while (<IN>) {
	chomp;
	my ($clusterId, $ctgNum, $cluster_len, $robustness, $cluster_str) = split /\s+/;
	$ClusterLen{$clusterId} = $cluster_len;
	my @t = split /;/, $cluster_str;
	foreach my $ctgId (@t) {
		$ctgId =~ s/[+-]//;
		$CtgCluster{$ctgId} = $clusterId;
	}
}
close IN;

#print Dumper \%CtgCluster;


##ptg000007l      ptg000018l      84733.00        head    tail
open IN,  $halfContig_contact_file || die "fail  $halfContig_contact_file";
while (<IN>) {
	next if(/^\#/);
	my ($ctg1,$ctg2,$contact) = split /\s+/;
	if (exists $CtgCluster{$ctg1} && exists $CtgCluster{$ctg2}) {
		my $cluster1 = $CtgCluster{$ctg1};
		my $cluster2 = $CtgCluster{$ctg2};
		if ($cluster1 lt $cluster2) {
			$ClusterContact{$cluster1}{$cluster2} += $contact;
		}else{
			$ClusterContact{$cluster2}{$cluster1} += $contact;
		}
	}
}	
close IN;

#print Dumper \%ClusterContact;


##print Dumper \%ClusterLen;

my $max_value = 0;
my @Output;
foreach my $cluster1 (sort keys %ClusterContact) {
	my $cluster1_p = $ClusterContact{$cluster1};
	foreach my $cluster2 (sort keys %$cluster1_p) {
		next if($cluster1 eq $cluster2);
		
		my $contact_normalized = $cluster1_p->{$cluster2} / ($ClusterLen{$cluster1} * $ClusterLen{$cluster2} / 1000000 / 1000000 ); 
		if ($max_value < $contact_normalized) {
			$max_value = $contact_normalized;
		}
		push  @Output, [$cluster1,$cluster2,$contact_normalized] ;
	}
}

print STDERR "MAX value : $max_value\n";

print "#cluster1\tcluster2\tdistance[0-1]\tcontact_normalized_scaled[max:1]\tcontact_normalized\n";
foreach my $p (@Output) {
	my ($cluster1,$cluster2,$contact_normalized) = @$p;
	my $contact_normalized_scaled = $contact_normalized / $max_value;
	my $distance = 1 - $contact_normalized_scaled;
	print "$cluster1\t$cluster2\t$distance\t$contact_normalized_scaled\t$contact_normalized\n";
}

#print Dumper \%ClusterContact;