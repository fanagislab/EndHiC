#!/usr/bin/perl

=head1 Name

summarize_endhic_results.pl  -- summarize and analyze a set of EndHiC results 

=head1 Description

Perform robustness analysis for all cluster types and find out stable cluster types

=head1 Version

  Author: Fan Wei, fanwei@caas.cn
  Version: 1.0,  Date: 2021/7/23
  Note:

=head1 Usage
  
  summarize_endhic_results.pl [options] <*.reciprocalMax.gfa.cluster.order.orient>
  --ctglenfile <str>  contig id and length file
  --ctglencut <int>   length cutoff for used large contigs, default=0
  --binsize <int>   bin size used in the hicpro result files, default=100000
  --help            output help information to screen  

=head1 Example

  perl summerize_endhic_results.pl  --binsize 100000  *.reciprocalMax.gfa.cluster.order.orient  > z.EndHiC.results.summary.and.analysist

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;


##get options from command line into variables and set default values
my $Contig_used_file;
my $Contig_len_cutoff;
my $Times;
my $BinSize;
my ($Verbose,$Help);
GetOptions(
	"ctglenfile:s"=>\$Contig_used_file,
	"ctglencut:i"=>\$Contig_len_cutoff,
	"binsize:i"=>\$BinSize,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$BinSize ||= 100000;
$Contig_len_cutoff ||= 0;
die `pod2text $0` if (@ARGV == 0 || $Help);


my $ref_cluster_file = shift;
my $hcluster_log_file = shift;

my %CtgChr;


##chr10   ptg000040l      ptg000003l      ptg000073l
open IN,  $ref_cluster_file || die "fail  $ref_cluster_file";
while (<IN>) {
	chomp;
	my @t = split /\s+/;
	my $chrId = shift @t;
	my $ctgnum = shift @t;
	foreach  my $ctgId (@t) {
		$CtgChr{$ctgId}  = $chrId;
	}
}
close IN;

#print Dumper \%CtgChr;


##Cluster_01      4       242687818       90      ptg000014l-;ptg000033l+;ptg000050l-;ptg000021l-
open IN,  $hcluster_log_file || die "fail  $hcluster_log_file";
while (<IN>) {
	if(!/^Cluster_/){
		print $_;
		next;
	}
	
	chomp;
	my @t = split /\s+/;
	my $chrId = shift @t;
	my $ctgnum = shift @t;
	my $ctgs_str = join(";",@t);
	my $chr_str;
	my %ChrCount;
	foreach  my $ctgId (@t) {
		my $chrId = $CtgChr{$ctgId};
		$chr_str .= "$chrId;";
		$ChrCount{$chrId} ++;
	}
	my $chr_num = keys %ChrCount;
	
	print "$chrId\t$ctgnum\t$ctgs_str\t$chr_num\t$chr_str\n";

}
close IN;

