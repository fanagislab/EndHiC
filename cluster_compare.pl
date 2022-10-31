#!/usr/bin/perl

=head1 Name

cluster_compare.pl  -- compare two cluster format files

=head1 Description 


=head1 Version

  Author: Fan Wei, fanwei@caas.cn
  Version: 1.0,  Date: 2021/7/23
  Note:

=head1 Usage
  
  cluster_compare.pl  <cluster_file1> <cluster_file2>

  --help            output help information to screen  

=head1 Example


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


my $cluster_file1 = shift;
my $cluster_file2 = shift;

my %CtgChr;
my %ChrCtg;


open IN,  $cluster_file1 || die "fail  $cluster_file1";
while (<IN>) {
	chomp;
	next if(/^\#/);
	my ($chrId, $ctgNum, $chr_len, $robustness, $ctgs) = split /\s+/;

	$ChrCtg{$chrId} = $ctgs;
	my @t = split /;/, $ctgs;

	foreach  my $ctgId (@t) {
		$ctgId =~ s/[+-]//;
		$CtgChr{$ctgId}  = $chrId;
	}
}
close IN;

#print Dumper \%CtgChr;

my %ChrClusters; 

##Cluster_01      4       242687818       90      ptg000014l-;ptg000033l+;ptg000050l-;ptg000021l-
open IN,  $cluster_file2 || die "fail  $cluster_file2";
while (<IN>) {
	my $line = $_;
	chomp;
	print $_;

	next if(/^\#/);
	my ($clusterId, $ctgNum, $cluster_len, $robustness, $cluster_str) = split /\s+/;
	
	#print STDERR $clusterId."\t$cluster_str\n";
	my @t = split /;/, $cluster_str;
	
	my %Chr;
	my $str;
	foreach my $ctgId (@t) {
		$ctgId =~ s/[+-]//;
		
		my $chrId = $CtgChr{$ctgId};
		#print STDERR "\t$ctgId\t$chrId\n";
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