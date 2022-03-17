#!/usr/bin/perl

=head1 Name

endhic_auto.pl  -- automatically run multiple rounds of endhic.pl

=head1 Description

Automatically run multiple rounds of endhic.pl, each round with parameters:

	A: " --minbinnum 5  --maxbinnum 25 "
	B: " --minbinnum 30 --maxbinnum 50 "
	C: " --minbinnum 55 --maxbinnum 75 "
	D: " --minbinnum 80 --maxbinnum 100 "
	E: " --minbinnum 105 --maxbinnum 125 "
	F: " --minbinnum 130 --maxbinnum 150 "
	G: " --minbinnum 155 --maxbinnum 175 "
	H: " --minbinnum 180 --maxbinnum 200 "
	.....................................
	.....................................
	.....................................

You can choose the best result from all rounds, in which the cluster number equals or approximates  the known chromosome number. Note that the assembly accuracy will decrease fastly as the running rounds growing. 

=head1 Version

  Author: Fan Wei, fanweiagis@126.com
  Version: 1.0,  Date: 2022/3/17

=head1 Usage

  --rounds <int>  number of running rounds of endhic.pl , range 1 - 8, default = 3
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

  perl ~/fanwei/code/AGIS/EndHiC/endhic_auto.pl --rounds 3  contigs_all.len combined_100000_abs.bed combined_100000.matrix combined_100000_iced.matrix


=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;


##get options from command line into variables and set default values
my ($RoundNum, $Verbose,$Help);
GetOptions(
	"rounds:i"=>\$RoundNum,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$RoundNum ||= 3;
die `pod2text $0` if (@ARGV == 0 || $Help);


my @Rounds = ("A", "B", "C", "D", "E", "F", "G", "H");
my @Parameters = (" --minbinnum 5  --maxbinnum 25 ", " --minbinnum 30 --maxbinnum 50 ", " --minbinnum 55 --maxbinnum 75 ", " --minbinnum 80 --maxbinnum 100 ", " --minbinnum 105 --maxbinnum 125 ", " --minbinnum 130 --maxbinnum 150 ", " --minbinnum 155 --maxbinnum 175 ", " --minbinnum 180 --maxbinnum 200 ");

die "\nError: exceed max allowed rounds\n\n" if($RoundNum > @Rounds);

my $contig_len_file = shift;
my $hicpro_bed_file = shift;
my $hicpro_raw_matrix_file = shift;
my $hicpro_iced_matrix_file = shift;


print STDERR "run round $Rounds[0]\n";
`perl $Bin/endhic.pl --outprefix Round_$Rounds[0]   $contig_len_file $hicpro_bed_file $hicpro_raw_matrix_file $hicpro_iced_matrix_file;  ln -s Round_$Rounds[0].04.summary_and_merging_results/z.EndHiC.$Rounds[0].results.summary.cluster ./`;

for (my $i = 1; $i < $RoundNum; $i ++) {
	
	my $current_mark = $Rounds[$i];
	print STDERR "run round $current_mark\n";
	my $last_mark = $Rounds[$i-1];
	`perl $Bin/cluster2bed.pl $hicpro_bed_file z.EndHiC.$last_mark.results.summary.cluster > cluster$last_mark\_100000_abs.bed  2> cluster$last_mark.id.len; perl $Bin/endhic.pl --outprefix Round_$Rounds[$i]  $Parameters[$i] --clustermark $current_mark cluster$last_mark.id.len  cluster$last_mark\_100000_abs.bed  $hicpro_raw_matrix_file $hicpro_iced_matrix_file; perl $Bin/cluster_merge.pl z.EndHiC.$last_mark.results.summary.cluster Round_$current_mark.04.summary_and_merging_results/z.EndHiC.$current_mark.results.summary.transit > Round_$current_mark.04.summary_and_merging_results/z.EndHiC.$current_mark.results.summary.cluster; ln -s Round_$current_mark.04.summary_and_merging_results/z.EndHiC.$current_mark.results.summary.cluster ./`;
	
}




####################################################
################### Sub Routines ###################
####################################################
