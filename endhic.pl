#!/usr/bin/perl

=head1 Name

endhic.pl  -- the standard pipeline of EndHic programs

=head1 Description

EndHiC calculates scaffolding results with various contig ends sizes, using various contact cutoff, for both raw and normalized contact data.

Then, EndHiC performs robustness analysis for all the cluster types in all the cluster results, and find the stable cluster types with dominant frequency.

Using default parameters:  "--binsize 100000  --minbinnum 5 --maxbinnum 25 --binnumstep 5",  represents using contig end sizes 0.5, 1, 1.5, 2, 2.5 Mb

Automatically using recommended sets of contact cutoffs: 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5 times of the turninging point. 

=head1 Version

  Author: Fan Wei, fanwei@caas.cn
  Version: 1.0,  Date: 2021/7/23
  Note:

=head1 Usage
  
  endhic.pl [options] <contig_id_len_file>  <hicpro_bed_file> <hicpro_raw_matrix_file> <hicpro_iced_matrix_file>
  --binsize <int>   bin size used in the hicpro result files, default=100000
  --minbinnum <int>    min bin number for contig end, head or tail end, default=5 
  --maxbinnum <int>      max bin number for contig end, head or tail end, default=25 
  --binnumstep <int>     step for the bin number increasing, default=5
  --clustermark <str>    mark character for cluster ID, default=A
  --outprefix <str>    prefix for the output subdirectories
  --help            output help information to screen  

=head1 Example
     
     perl endhic.pl  contigs.id.len formal_100000_abs.bed formal_100000.matrix formal_1000000_iced.matrix


=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;


##get options from command line into variables and set default values
my $BinSize;
my $MinBinNum;
my $MaxBinNum;
my $BinNumStep;
my $ClusterMark;
my $expected_chromosome_number;  
my $minimum_chromosome_length;
my $TimesStart;
my $TimesStop;
my $TimesStep;
my $OutPrefix;
my ($Verbose,$Help);
GetOptions(
	"binsize:i"=>\$BinSize,
	"minbinnum:i"=>\$MinBinNum,
	"maxbinnum:i"=>\$MaxBinNum,
	"binnumstep:i"=>\$BinNumStep,
	"clustermark:s"=>\$ClusterMark,
	"chrnum:i"=>\$expected_chromosome_number,  ##hidden parameter, not for users
	"minchrlen:i"=>\$minimum_chromosome_length, ##hidden parameter, not for users
	"mintimes:f"=>\$TimesStart, ##hidden parameter, not for users
	"maxtimes:f"=>\$TimesStop, ##hidden parameter, not for users
	"timesstep:f"=>\$TimesStep, ##hidden parameter, not for users
	"outprefix:s"=>\$OutPrefix,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$BinSize ||= 100000;
$MinBinNum ||= 5;
$MaxBinNum ||= 25;
$BinNumStep ||= 5;
$ClusterMark ||= "A";
$expected_chromosome_number ||= 0;
$minimum_chromosome_length ||= 10000000;
$TimesStart ||= 1.5;
$TimesStop  ||= 4.5;
$TimesStep  ||= 0.5;
$OutPrefix ||= "Round_A";
die `pod2text $0` if (@ARGV == 0 || $Help);



my $contig_len_file = shift;
my $hicpro_bed_file = shift;
my $hicpro_raw_matrix_file = shift;
my $hicpro_iced_matrix_file = shift;


my $work_shell_file = "endhic.$OutPrefix.sh";
open ENDHIC, ">$work_shell_file" || die "fail open $work_shell_file\n";


##EndHic for fixed-size contig ends
my $job_num = 0;
for (my $BinNum = $MinBinNum; $BinNum <= $MaxBinNum; $BinNum += $BinNumStep) {
	print ENDHIC "perl $Bin/endhic_ctgEnd_pipeline.pl --binsize $BinSize --binnum $BinNum --mintimes $TimesStart  --maxtimes $TimesStop  --timesstep  $TimesStep  --chrnum $expected_chromosome_number --minchrlen $minimum_chromosome_length  $contig_len_file $hicpro_bed_file $hicpro_raw_matrix_file\n";
	print ENDHIC "perl $Bin/endhic_ctgEnd_pipeline.pl --binsize $BinSize --binnum $BinNum --mintimes $TimesStart  --maxtimes $TimesStop  --timesstep  $TimesStep  --chrnum $expected_chromosome_number --minchrlen $minimum_chromosome_length  $contig_len_file $hicpro_bed_file $hicpro_iced_matrix_file\n";
	$job_num += 2;
}

close ENDHIC;

##run the work shell in parallel
`perl $Bin/multi-process.pl -cpu $job_num  $work_shell_file`;

open ENDHIC, ">>$work_shell_file" || die "fail open $work_shell_file\n";

print ENDHIC "\n";



##summarize the major EndHic results
my $CtgLenCut = $BinSize * $MaxBinNum * 2; 
my $final_file_end = ($ClusterMark eq "A") ? "cluster" : "transit";

print ENDHIC "perl $Bin/summarize_endhic_results.pl --clustermark $ClusterMark  --binsize $BinSize  --ctglencut $CtgLenCut  --ctglenfile $contig_len_file  $hicpro_raw_matrix_file.*.gfa.cluster  $hicpro_iced_matrix_file.*.gfa.cluster  > z.EndHiC.$ClusterMark.results.summary  2> z.EndHiC.$ClusterMark.results.summary.$final_file_end\n";
             `perl $Bin/summarize_endhic_results.pl --clustermark $ClusterMark  --binsize $BinSize  --ctglencut $CtgLenCut  --ctglenfile $contig_len_file  $hicpro_raw_matrix_file.*.gfa.cluster  $hicpro_iced_matrix_file.*.gfa.cluster  > z.EndHiC.$ClusterMark.results.summary  2> z.EndHiC.$ClusterMark.results.summary.$final_file_end`;

`perl $Bin/cluster_to_GFA.pl $contig_len_file z.EndHiC.$ClusterMark.results.summary.$final_file_end > z.EndHiC.$ClusterMark.results.summary.$final_file_end.GFA`;

`perl $Bin/generate_GFAv1.2.pl z.EndHiC.$ClusterMark.results.summary.$final_file_end.GFA > z.EndHiC.$ClusterMark.results.summary.$final_file_end.GFA.v1.2.GFA`;

print ENDHIC "perl $Bin/summarize_endhic_results.pl --clustermark $ClusterMark  --binsize $BinSize  --ctglencut $CtgLenCut  --ctglenfile $contig_len_file  $hicpro_raw_matrix_file.*.gfa.cluster    > z.EndHiC.$ClusterMark.results.raw.summary  2> z.EndHiC.$ClusterMark.results.raw.summary.$final_file_end\n";
             `perl $Bin/summarize_endhic_results.pl --clustermark $ClusterMark  --binsize $BinSize  --ctglencut $CtgLenCut  --ctglenfile $contig_len_file  $hicpro_raw_matrix_file.*.gfa.cluster    > z.EndHiC.$ClusterMark.results.raw.summary  2> z.EndHiC.$ClusterMark.results.raw.summary.$final_file_end`;

`perl $Bin/cluster_to_GFA.pl $contig_len_file z.EndHiC.$ClusterMark.results.raw.summary.$final_file_end > z.EndHiC.$ClusterMark.results.raw.summary.$final_file_end.GFA`;

`perl $Bin/generate_GFAv1.2.pl z.EndHiC.$ClusterMark.results.raw.summary.$final_file_end.GFA > z.EndHiC.$ClusterMark.results.raw.summary.$final_file_end.GFA.v1.2.GFA`;


print ENDHIC "perl $Bin/summarize_endhic_results.pl --clustermark $ClusterMark  --binsize $BinSize  --ctglencut $CtgLenCut  --ctglenfile $contig_len_file  $hicpro_iced_matrix_file.*.gfa.cluster  > z.EndHiC.$ClusterMark.results.iced.summary  2> z.EndHiC.$ClusterMark.results.iced.summary.$final_file_end\n";
             `perl $Bin/summarize_endhic_results.pl --clustermark $ClusterMark  --binsize $BinSize  --ctglencut $CtgLenCut  --ctglenfile $contig_len_file  $hicpro_iced_matrix_file.*.gfa.cluster  > z.EndHiC.$ClusterMark.results.iced.summary  2> z.EndHiC.$ClusterMark.results.iced.summary.$final_file_end`;

`perl $Bin/cluster_to_GFA.pl $contig_len_file z.EndHiC.$ClusterMark.results.iced.summary.$final_file_end > z.EndHiC.$ClusterMark.results.iced.summary.$final_file_end.GFA`;

`perl $Bin/generate_GFAv1.2.pl z.EndHiC.$ClusterMark.results.iced.summary.$final_file_end.GFA > z.EndHiC.$ClusterMark.results.iced.summary.$final_file_end.GFA.v1.2.GFA`;

if ($expected_chromosome_number > 0) {
	print ENDHIC "perl $Bin/summarize_endhic_results.pl --clustermark $ClusterMark  --binsize $BinSize  --ctglencut $CtgLenCut  --ctglenfile $contig_len_file   $hicpro_raw_matrix_file.*.hcluster.need   $hicpro_iced_matrix_file.*.hcluster.need   > z.EndHiC.$ClusterMark.results.verify.summary  2> z.EndHiC.$ClusterMark.results.verify.summary.$final_file_end\n";
				 `perl $Bin/summarize_endhic_results.pl --clustermark $ClusterMark  --binsize $BinSize  --ctglencut $CtgLenCut  --ctglenfile $contig_len_file   $hicpro_raw_matrix_file.*.hcluster.need   $hicpro_iced_matrix_file.*.hcluster.need   > z.EndHiC.$ClusterMark.results.verify.summary  2> z.EndHiC.$ClusterMark.results.verify.summary.$final_file_end`;

}



##re-organize the output files, move into subdirectories
print ENDHIC "mkdir $OutPrefix.01.contig_end_contact_results   $OutPrefix.02.GFA_contig_graph_results  $OutPrefix.03.cluster_order_orient_results  $OutPrefix.04.summary_and_merging_results\n";
             `mkdir $OutPrefix.01.contig_end_contact_results   $OutPrefix.02.GFA_contig_graph_results  $OutPrefix.03.cluster_order_orient_results  $OutPrefix.04.summary_and_merging_results `;

print ENDHIC "mv z.EndHiC.$ClusterMark.*  $hicpro_raw_matrix_file.*.summary* $hicpro_iced_matrix_file.*.summary*   $OutPrefix.04.summary_and_merging_results\n";
             `mv z.EndHiC.$ClusterMark.*  $hicpro_raw_matrix_file.*.summary* $hicpro_iced_matrix_file.*.summary*   $OutPrefix.04.summary_and_merging_results`;

print ENDHIC "mv *.CtgContact *.adjustTransform  *.turningPoint  $OutPrefix.01.contig_end_contact_results/\n";
             `mv *.CtgContact *.adjustTransform  *.turningPoint  $OutPrefix.01.contig_end_contact_results/`;

print ENDHIC "mv *.gfa  $OutPrefix.02.GFA_contig_graph_results/\n";
             `mv *.gfa  $OutPrefix.02.GFA_contig_graph_results/`;

print ENDHIC "mv *.gfa.cluster  *.gfa.topology  $OutPrefix.03.cluster_order_orient_results\n";
             `mv *.gfa.cluster  *.gfa.topology  $OutPrefix.03.cluster_order_orient_results`;


if ($expected_chromosome_number > 0) {
	print ENDHIC "mkdir 05.hierarchical_cluster_verify;  mv  *.CtgContact.distance  *.distance.hcluster  *.hcluster.need  05.hierarchical_cluster_verify/\n";
				 `mkdir 05.hierarchical_cluster_verify;  mv  *.CtgContact.distance  *.distance.hcluster  *.hcluster.need  05.hierarchical_cluster_verify/`;
}


close ENDHIC;



####################################################
################### Sub Routines ###################
####################################################
