#!/usr/bin/perl

=head1 Name

endhic_iterate.pl  -- automatically run multiple iterative rounds of endhic.pl

=head1 Description

Automatically run multiple rounds (max:26) of endhic.pl, each round with increasing contig end sizes. 

Below is an example when the parameter "--binnumstep" is set to 5:
	
	A: " --minbinnum 5  --maxbinnum 25 --binnumstep 5 "

	B: " --minbinnum 30 --maxbinnum 50 --binnumstep 5 "

	C: " --minbinnum 55 --maxbinnum 75 --binnumstep 5 "

	D: " --minbinnum 80 --maxbinnum 100 --binnumstep 5 "

	E: " --minbinnum 105 --maxbinnum 125 --binnumstep 5 "

	.....................................
	.....................................
	.....................................

Below is an example when the parameter "--binnumstep" is set to 1:
	
	A: " --minbinnum 5  --maxbinnum 9 --binnumstep 1 "
	
	B: " --minbinnum 10 --maxbinnum 14 --binnumstep 1 "
	
	C: " --minbinnum 15 --maxbinnum 19 --binnumstep 1 "
	
	D: " --minbinnum 20 --maxbinnum 24 --binnumstep 1 "
	
	E: " --minbinnum 25 --maxbinnum 29 --binnumstep 1 "
	
	.....................................
	.....................................
	.....................................

For the meaning of these parameters, please see the command-line help of endhic.pl

You can choose the best result from all rounds, in which the cluster number equals or approximates  the known chromosome number. Note that the assembly accuracy will decrease as the running rounds growing. 

=head1 Version

  Author: Fan Wei, fanweiagis@126.com
  Version: 1.0,  Date: 2022/3/17

=head1 Usage

  --rounds <int>  number of running rounds of endhic.pl , range 1 - 26, default = 3
  --binsize <int>   bin size used in the hicpro result files, default = 100000
  --binnumstep <int> the binnumstep parameter of endhic.pl, range 1 - 10, default = 5
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

  (1) For long-contiuous contig assembly, use default parameters
  endhic_iterate.pl --rounds 3 --binnumstep 5 contigs_all.len combined_100000_abs.bed combined_100000.matrix combined_100000_iced.matrix
  
  (2) For relatively shorter contigs, run more rounds with smaller increasing of contig end sizes 
  endhic_iterate.pl --rounds 5 --binnumstep 3 contigs_all.len combined_100000_abs.bed combined_100000.matrix combined_100000_iced.matrix


=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;


##get options from command line into variables and set default values
my ($RoundNum, $Verbose,$BinSize, $binnumstep, $Help);
GetOptions(
	"rounds:i"=>\$RoundNum,
	"binsize:i"=>\$BinSize,
	"binnumstep:i"=>\$binnumstep,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$RoundNum ||= 3;
$BinSize ||= 100000;
$binnumstep ||= 5;
die `pod2text $0` if (@ARGV == 0 || $Help);

my @Rounds;
my @Parameters;
#my @Rounds = ("A", "B", "C", "D", "E", "F", "G", "H");
#my @Parameters = (" --minbinnum 5  --maxbinnum 25 ", " --minbinnum 30 --maxbinnum 50 ", " --minbinnum 55 --maxbinnum 75 ", " --minbinnum 80 --maxbinnum 100 ", " --minbinnum 105 --maxbinnum 125 ", " --minbinnum 130 --maxbinnum 150 ", " --minbinnum 155 --maxbinnum 175 ", " --minbinnum 180 --maxbinnum 200 ");


##generate the @Rounds and @Parameters
my $MaxRounds = 26;  ##Rounds A-Z
my $character = "A";
for (my $i = 0; $i < $MaxRounds; $i++) {
	my $minbinnum = 5 + $i * 5 * $binnumstep;
	my $maxbinnum = 5 + ($i + 1) * 5 * $binnumstep - $binnumstep;
	
	push @Rounds, $character;
	push @Parameters, " --binsize $BinSize --minbinnum $minbinnum  --maxbinnum $maxbinnum  --binnumstep $binnumstep ";
	$character ++;
}


die "\nError: exceed max allowed rounds\n\n" if($RoundNum > @Rounds);

my $contig_len_file = shift;
my $hicpro_bed_file = shift;
my $hicpro_raw_matrix_file = shift;
my $hicpro_iced_matrix_file = shift;

open OUT, ">endhic_iterate.sh" || die "fail";

print STDERR "run round $Rounds[0]\n";
`perl $Bin/endhic.pl --outprefix Round_$Rounds[0] $Parameters[0] --clustermark $Rounds[0]   $contig_len_file $hicpro_bed_file $hicpro_raw_matrix_file $hicpro_iced_matrix_file;  ln -s Round_$Rounds[0].04.summary_and_merging_results/z.EndHiC.$Rounds[0].results.summary.cluster* ./`;

print OUT "perl $Bin/endhic.pl --outprefix Round_$Rounds[0] $Parameters[0] --clustermark $Rounds[0]   $contig_len_file $hicpro_bed_file $hicpro_raw_matrix_file $hicpro_iced_matrix_file;  ln -s Round_$Rounds[0].04.summary_and_merging_results/z.EndHiC.$Rounds[0].results.summary.cluster* ./\n\n";

for (my $i = 1; $i < $RoundNum; $i ++) {
	
	my $current_mark = $Rounds[$i];
	print STDERR "run round $current_mark\n";
	my $last_mark = $Rounds[$i-1];
	`perl $Bin/cluster2bed.pl $hicpro_bed_file z.EndHiC.$last_mark.results.summary.cluster > cluster$last_mark\_100000_abs.bed  2> cluster$last_mark.id.len; perl $Bin/endhic.pl --outprefix Round_$Rounds[$i]  $Parameters[$i] --clustermark $current_mark cluster$last_mark.id.len  cluster$last_mark\_100000_abs.bed  $hicpro_raw_matrix_file $hicpro_iced_matrix_file; perl $Bin/cluster_merge.pl z.EndHiC.$last_mark.results.summary.cluster Round_$current_mark.04.summary_and_merging_results/z.EndHiC.$current_mark.results.summary.transit > Round_$current_mark.04.summary_and_merging_results/z.EndHiC.$current_mark.results.summary.cluster; perl $Bin/cluster_to_GFA.pl $contig_len_file Round_$current_mark.04.summary_and_merging_results/z.EndHiC.$current_mark.results.summary.cluster > Round_$current_mark.04.summary_and_merging_results/z.EndHiC.$current_mark.results.summary.cluster.GFA;   perl $Bin/generate_GFAv1.2.pl Round_$current_mark.04.summary_and_merging_results/z.EndHiC.$current_mark.results.summary.cluster.GFA > Round_$current_mark.04.summary_and_merging_results/z.EndHiC.$current_mark.results.summary.cluster.GFA.v1.2.GFA;  ln -s Round_$current_mark.04.summary_and_merging_results/z.EndHiC.$current_mark.results.summary.cluster* ./`;

	print OUT "perl $Bin/cluster2bed.pl $hicpro_bed_file z.EndHiC.$last_mark.results.summary.cluster > cluster$last_mark\_100000_abs.bed  2> cluster$last_mark.id.len;\n perl $Bin/endhic.pl --outprefix Round_$Rounds[$i]  $Parameters[$i] --clustermark $current_mark cluster$last_mark.id.len  cluster$last_mark\_100000_abs.bed  $hicpro_raw_matrix_file $hicpro_iced_matrix_file;\n perl $Bin/cluster_merge.pl z.EndHiC.$last_mark.results.summary.cluster Round_$current_mark.04.summary_and_merging_results/z.EndHiC.$current_mark.results.summary.transit > Round_$current_mark.04.summary_and_merging_results/z.EndHiC.$current_mark.results.summary.cluster;\n perl $Bin/cluster_to_GFA.pl $contig_len_file Round_$current_mark.04.summary_and_merging_results/z.EndHiC.$current_mark.results.summary.cluster > Round_$current_mark.04.summary_and_merging_results/z.EndHiC.$current_mark.results.summary.cluster.GFA;\n  perl $Bin/generate_GFAv1.2.pl Round_$current_mark.04.summary_and_merging_results/z.EndHiC.$current_mark.results.summary.cluster.GFA > Round_$current_mark.04.summary_and_merging_results/z.EndHiC.$current_mark.results.summary.cluster.GFA.v1.2.GFA; \n ln -s Round_$current_mark.04.summary_and_merging_results/z.EndHiC.$current_mark.results.summary.cluster* ./\n\n";
	
}

close OUT;

####################################################
################### Sub Routines ###################
####################################################
