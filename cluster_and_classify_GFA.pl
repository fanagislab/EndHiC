#!/usr/bin/perl

=head1 Name

cluster_and_classify.pl  --  cluster the GFA using hcluter algorithm and distinguish linear and circular topology

=head1 Description

This program invoke "hcluster_from_GFA.pl" to perform clustering function

The major code of this program focus on distinguishing the linear and circular topology


=head1 Version

  Author: Fan Wei, fanw@caas.cn
  Version: 1.0,  Date: 2021/8/22
  Note:

=head1 Usage

  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

  cluster_and_classify_GFA.pl formal_1000000_iced.matrix.revised.CtgContactFilter.overCutoff.5.0.reciprocalMax.gfa > formal_1000000_iced.matrix.revised.CtgContactFilter.overCutoff.5.0.reciprocalMax.gfa.cluster

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

my $gfa_file = shift;

my $hcluster_res_file = "$gfa_file.hcluster";

`perl $Bin/hcluster_from_GFA.pl $gfa_file > $hcluster_res_file`;


my %Links;


open IN, $gfa_file || die "fail open $gfa_file";
while (<IN>) {
	
	my $line = $_;

	if (/^L/) {
		#L       ptg000002l      +       ptg000026l      +       0M      L1:i:0  RM:i:1
		my @t = split /\s+/;
		my $type = $t[0];
		my $ctg1 = $t[1];
		my $ctg1_strand = $t[2];
		my $ctg2 = $t[3];
		my $ctg2_strand = $t[4];
		my $match = $t[5];
		my $L1 = $t[6];
		my $RM = $t[7];
		my $contact = $t[8];

		push @{$Links{$ctg1}}, $ctg2;
		push @{$Links{$ctg2}}, $ctg1;

		
	}
}
close IN;


my @Output;
open IN, $hcluster_res_file || die "fail open $hcluster_res_file";
while (<IN>) {
	my @t = split /\s+/;
	
	my $cluster_id = $t[0];
	my $ctgs_num = $t[1];
	my @ctgs;
	for (my $i=2; $i<@t; $i++) {
		push @ctgs, $t[$i];
	}

	my $LinkOne_count = 0;
	my $LinkTwo_count = 0;
	my $LinkMore_count = 0;

	foreach my $ctg (@ctgs) {
		my $ctg_p = $Links{$ctg};
		my $ctg_links = @$ctg_p;
		$LinkOne_count ++ if($ctg_links == 1);
		$LinkTwo_count ++ if($ctg_links == 2);
		$LinkMore_count ++ if($ctg_links > 2);
	}
	#print "$cluster_id\t$LinkOne_count\t$LinkTwo_count\n";

	my $Class;
	if ($LinkOne_count == 2 && $LinkTwo_count == $ctgs_num - 2) {
		$Class = "Linear";
	}
	if ($LinkOne_count == 0 && $LinkTwo_count == $ctgs_num ) {
		$Class = "Circular";
	}
	if ($LinkMore_count > 0) {
		$Class = "Branch";
	}
	#print "$cluster_id\t$LinkOne_count\t$LinkTwo_count\t$Class\n";
	push @t, "[Type:$Class]";
	my $line = join("\t",@t);
	#print $line."\n";
	push @Output, $line;

}
close IN;


#output the results
foreach my $line (@Output) {
	print  $line."\n";
}

`rm $hcluster_res_file`;



