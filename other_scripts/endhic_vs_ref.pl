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


my $ref_contig_pos_file = shift;
my $endhic_result_file = shift;

my %CtgChr;
my %ChrCtg;


open OUT, ">$ref_contig_pos_file.contigs" || die "fail $ref_contig_pos_file.contigs";
open IN,  $ref_contig_pos_file || die "fail  $ref_contig_pos_file";
while (<IN>) {
	chomp;
	my @t = split /\s+/;
	my $chrId = shift @t;
	my $ctgs = join("  ",@t);
	$ChrCtg{$chrId} = $ctgs;
	foreach  my $ctgId (@t) {
		$CtgChr{$ctgId}  = $chrId;
		print OUT "$ctgId\t$chrId\n";
	}
}
close IN;
close OUT;

#print Dumper \%CtgChr;

my %ChrClusters; 

##Cluster_01      4       242687818       90      ptg000014l-;ptg000033l+;ptg000050l-;ptg000021l-
open IN,  $endhic_result_file || die "fail  $endhic_result_file";
while (<IN>) {
	my $line = $_;
	chomp;
	print $_;
	my ($clusterId, $ctgNum, $cluster_len, $robustness, $cluster_str) = split /\s+/;
	
	print STDERR $clusterId."\t$cluster_str\n";
	my @t = split /;/, $cluster_str;
	
	my %Chr;
	my $str;
	foreach my $ctgId (@t) {
		$ctgId =~ s/[+-]//;
		
		my $chrId = $CtgChr{$ctgId};
		print STDERR "\t$ctgId\t$chrId\n";
		$str .= "$chrId;";
		$Chr{$chrId} ++;
	}
	print "\t".$str."\t";

	my @Chr = keys %Chr;
	if (@Chr > 1) {
		print STDERR "EndHiC scaffolding error at $clusterId\n";
	}elsif(@Chr == 1 ){
		print $Chr[0]."\t";
		push @{$ChrClusters{$Chr[0]}}, $line;
	}else{
		print STDERR "EndHiC scaffolding error at $clusterId, no this contig\n";
	}
	print "\n";


}
close IN;


print "\n--------------------------------------------------\n";

foreach my $chr (sort keys %ChrClusters) {
	my $chr_p = $ChrClusters{$chr};
	print "\n".$chr.":  $ChrCtg{$chr}\n";
	foreach my $line (@$chr_p) {
		print $line;
	}
}