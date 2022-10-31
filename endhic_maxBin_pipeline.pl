#!/usr/bin/perl

=head1 Name

endhic_maxBin_pipeline.pl  -- the pipeline of endhic using the max bin algorithm

=head1 Description

Invoke these programs sequentially: 
	(1) revise_iced_matrix.pl                check the iced matrix and recover those with raw contact lower than a given cutoff
	(2) ctgContact_from_maxBinContacts.pl   calculate the max bin contact for each contig pairs; 
	(3) filter_CtgContact.pl           filter the contig contacts by bin distance to head or tail ends, and sort 
	(4) scaffold_by_trueCtgContact.pl        scaffolding by link the contigs with true HiC contacts, output GFA file;
	(5) turningpoint_by_lineartransform.pl  get the turning point by adjusting the contacts data and linear transforming
	(6) cluster_and_classify_GFA.pl         cluster from GFA and distinguish linear and circular topology
	(7) order_and_orient_GFA.pl             get order and orientation for each cluster from GFA
	(8) summerize_endhic_results.pl         summarize a set of Endhic results

=head1 Version

  Author: Fan Wei, fanwei@caas.cn
  Version: 1.0,  Date: 2021/7/23
  Note:

=head1 Usage
  
  endhic_maxBin_pipeline.pl [options] <contig_id_len_file>  <hicpro_bed_file> <hicpro_matrix_file> 
  --binsize <int>   bin size for contig ends (head and tail)
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

  perl endhic_maxBin_pipeline.pl --binsize 1000000  contigs.id.len formal_1000000_abs.bed formal_1000000.matrix
  perl endhic_maxBin_pipeline.pl --binsize 1000000  contigs.id.len formal_1000000_abs.bed formal_1000000_iced.matrix

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;


##get options from command line into variables and set default values
my $BinSize;
my $Times;
my ($Verbose,$Help);
GetOptions(
	"binsize:i"=>\$BinSize,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$BinSize ||= 1000000;
die `pod2text $0` if (@ARGV == 0 || $Help);

my $contig_len_file = shift;
my $hicpro_bed_file = shift;
my $hicpro_matrix_file = shift;

my $endhic_result_files;

print STDERR "Using Bin Size:  $BinSize\n";
print STDERR "Times of the turning point to be used as reliable contacts cutoff:  $Times\n";

my $CtgSize_cutoff = $BinSize * 2;
my $rawcontacts_cutoff = $BinSize / 10000;

print STDERR "Contigs with length shorter than this cutoff were filtered in the first round of clustering:  $CtgSize_cutoff\n";
print STDERR "The raw contacts signal less than this cutoff will not be normalized in the iced matrix file: $rawcontacts_cutoff\n";

my $DistToEnd_cutoff;
if ($BinSize >= 1500000) {
	$DistToEnd_cutoff = 2;
}elsif ($BinSize >= 1000000) {
	$DistToEnd_cutoff = 3;
}elsif ($BinSize >= 500000) {
	$DistToEnd_cutoff = 6;
}elsif ($BinSize >= 200000) {
	$DistToEnd_cutoff = 9;
}elsif ($BinSize >= 100000) {
	$DistToEnd_cutoff = 15;
}else{
	$DistToEnd_cutoff = 20;
}
print STDERR "The max bin whose distance is larger than this cutoff to the head or tail bin will be filtered:  $DistToEnd_cutoff\n";

my $work_shell_file = "endhic.$BinSize.sh";

open ENDHIC, ">$work_shell_file" || die "fail open $work_shell_file\n";


########################################################################################################


print ENDHIC "perl $Bin/ctgContact_from_maxBinContacts.pl $contig_len_file $hicpro_bed_file $hicpro_matrix_file > $hicpro_matrix_file.CtgContact\n";
`perl $Bin/ctgContact_from_maxBinContacts.pl $contig_len_file $hicpro_bed_file $hicpro_matrix_file > $hicpro_matrix_file.CtgContact`;


print ENDHIC "$Bin/turningpoint_by_lineartransform.pl $hicpro_matrix_file.CtgContact > $hicpro_matrix_file.CtgContact.adjustTransform 2> $hicpro_matrix_file.CtgContact.turningPoint\n";
`perl $Bin/turningpoint_by_lineartransform.pl $hicpro_matrix_file.CtgContact > $hicpro_matrix_file.CtgContact.adjustTransform 2> $hicpro_matrix_file.CtgContact.turningPoint`;

my $raw_turning_point = `cat $hicpro_matrix_file.CtgContact.turningPoint`;
$raw_turning_point = $1  if($raw_turning_point =~ /Turning point:\s+(\S+)/);

print STDERR "All turning_point: $raw_turning_point\n";

my $raw_contact_cutoff;

for ($Times = 1.0; $Times <= 5.0; $Times += 0.5) {
	$Times = sprintf("%.1f",$Times);
	$raw_contact_cutoff = $Times * $raw_turning_point;

	print ENDHIC "perl $Bin/scaffold_by_trueCtgContact.pl   --contacts $raw_contact_cutoff $contig_len_file  $hicpro_matrix_file.CtgContact > $hicpro_matrix_file.CtgContact.overCutoff.$Times.gfa \n";
	`perl $Bin/scaffold_by_trueCtgContact.pl   --contacts $raw_contact_cutoff $contig_len_file  $hicpro_matrix_file.CtgContact > $hicpro_matrix_file.CtgContact.overCutoff.$Times.gfa `;

	print ENDHIC "perl $Bin/scaffold_by_trueCtgContact.pl   --reciprocalmax  $contig_len_file  $hicpro_matrix_file.CtgContact > $hicpro_matrix_file.CtgContact.reciprocalMax.gfa \n";
	`perl $Bin/scaffold_by_trueCtgContact.pl   --reciprocalmax  $contig_len_file  $hicpro_matrix_file.CtgContact > $hicpro_matrix_file.CtgContact.reciprocalMax.gfa `;

	print ENDHIC "perl $Bin/scaffold_by_trueCtgContact.pl   --contacts $raw_contact_cutoff --reciprocalmax  $contig_len_file  $hicpro_matrix_file.CtgContact > $hicpro_matrix_file.CtgContact.overCutoff.$Times.reciprocalMax.gfa \n";
	`perl $Bin/scaffold_by_trueCtgContact.pl   --contacts $raw_contact_cutoff --reciprocalmax  $contig_len_file  $hicpro_matrix_file.CtgContact > $hicpro_matrix_file.CtgContact.overCutoff.$Times.reciprocalMax.gfa `;

	print ENDHIC "perl $Bin/cluster_and_classify_GFA.pl $hicpro_matrix_file.CtgContact.overCutoff.$Times.reciprocalMax.gfa > $hicpro_matrix_file.CtgContact.overCutoff.$Times.reciprocalMax.gfa.topology\n";
	`perl $Bin/cluster_and_classify_GFA.pl $hicpro_matrix_file.CtgContact.overCutoff.$Times.reciprocalMax.gfa > $hicpro_matrix_file.CtgContact.overCutoff.$Times.reciprocalMax.gfa.topology`;

	print ENDHIC "perl $Bin/order_and_orient_GFA.pl --size $CtgSize_cutoff $hicpro_matrix_file.CtgContact.overCutoff.$Times.reciprocalMax.gfa  $hicpro_matrix_file.CtgContact.overCutoff.$Times.reciprocalMax.gfa.topology > $hicpro_matrix_file.CtgContact.overCutoff.$Times.reciprocalMax.gfa.cluster\n";
	`perl $Bin/order_and_orient_GFA.pl --size $CtgSize_cutoff $hicpro_matrix_file.CtgContact.overCutoff.$Times.reciprocalMax.gfa  $hicpro_matrix_file.CtgContact.overCutoff.$Times.reciprocalMax.gfa.topology > $hicpro_matrix_file.CtgContact.overCutoff.$Times.reciprocalMax.gfa.cluster`;

	$endhic_result_files .= "  $hicpro_matrix_file.CtgContact.overCutoff.$Times.reciprocalMax.gfa.cluster";
}

##################################################################

print ENDHIC "perl $Bin/filter_CtgContact.pl --disttoend  $DistToEnd_cutoff $contig_len_file  $hicpro_matrix_file.CtgContact > $hicpro_matrix_file.CtgContactFilter \n";
`perl $Bin/filter_CtgContact.pl  --disttoend  $DistToEnd_cutoff  $contig_len_file  $hicpro_matrix_file.CtgContact > $hicpro_matrix_file.CtgContactFilter `;

print ENDHIC "$Bin/turningpoint_by_lineartransform.pl $hicpro_matrix_file.CtgContactFilter > $hicpro_matrix_file.CtgContactFilter.adjustTransform 2> $hicpro_matrix_file.CtgContactFilter.turningPoint\n";
`perl $Bin/turningpoint_by_lineartransform.pl $hicpro_matrix_file.CtgContactFilter > $hicpro_matrix_file.CtgContactFilter.adjustTransform 2> $hicpro_matrix_file.CtgContactFilter.turningPoint`;

$raw_turning_point = `cat $hicpro_matrix_file.CtgContactFilter.turningPoint`;
$raw_turning_point = $1  if($raw_turning_point =~ /Turning point:\s+(\S+)/);
print STDERR "Filtered turning_point: $raw_turning_point\n";

for ($Times = 1.0; $Times <= 5.0; $Times += 0.5) {
	$Times = sprintf("%.1f",$Times);
	$raw_contact_cutoff = $Times * $raw_turning_point;

	print ENDHIC "perl $Bin/scaffold_by_trueCtgContact.pl   --contacts $raw_contact_cutoff $contig_len_file  $hicpro_matrix_file.CtgContactFilter > $hicpro_matrix_file.CtgContactFilter.overCutoff.$Times.gfa \n";
	`perl $Bin/scaffold_by_trueCtgContact.pl   --contacts $raw_contact_cutoff $contig_len_file  $hicpro_matrix_file.CtgContactFilter > $hicpro_matrix_file.CtgContactFilter.overCutoff.$Times.gfa `;

	print ENDHIC "perl $Bin/scaffold_by_trueCtgContact.pl   --reciprocalmax  $contig_len_file  $hicpro_matrix_file.CtgContactFilter > $hicpro_matrix_file.CtgContactFilter.reciprocalMax.gfa \n";
	`perl $Bin/scaffold_by_trueCtgContact.pl   --reciprocalmax  $contig_len_file  $hicpro_matrix_file.CtgContactFilter > $hicpro_matrix_file.CtgContactFilter.reciprocalMax.gfa `;

	print ENDHIC "perl $Bin/scaffold_by_trueCtgContact.pl   --contacts $raw_contact_cutoff --reciprocalmax  $contig_len_file  $hicpro_matrix_file.CtgContactFilter > $hicpro_matrix_file.CtgContactFilter.overCutoff.$Times.reciprocalMax.gfa \n";
	`perl $Bin/scaffold_by_trueCtgContact.pl   --contacts $raw_contact_cutoff --reciprocalmax  $contig_len_file  $hicpro_matrix_file.CtgContactFilter > $hicpro_matrix_file.CtgContactFilter.overCutoff.$Times.reciprocalMax.gfa `;

	print ENDHIC "perl $Bin/cluster_and_classify_GFA.pl $hicpro_matrix_file.CtgContactFilter.overCutoff.$Times.reciprocalMax.gfa > $hicpro_matrix_file.CtgContactFilter.overCutoff.$Times.reciprocalMax.gfa.topology\n";
	`perl $Bin/cluster_and_classify_GFA.pl $hicpro_matrix_file.CtgContactFilter.overCutoff.$Times.reciprocalMax.gfa > $hicpro_matrix_file.CtgContactFilter.overCutoff.$Times.reciprocalMax.gfa.topology`;

	print ENDHIC "perl $Bin/order_and_orient_GFA.pl --size $CtgSize_cutoff $hicpro_matrix_file.CtgContactFilter.overCutoff.$Times.reciprocalMax.gfa  $hicpro_matrix_file.CtgContactFilter.overCutoff.$Times.reciprocalMax.gfa.topology > $hicpro_matrix_file.CtgContactFilter.overCutoff.$Times.reciprocalMax.gfa.cluster\n";
	`perl $Bin/order_and_orient_GFA.pl --size $CtgSize_cutoff $hicpro_matrix_file.CtgContactFilter.overCutoff.$Times.reciprocalMax.gfa  $hicpro_matrix_file.CtgContactFilter.overCutoff.$Times.reciprocalMax.gfa.topology > $hicpro_matrix_file.CtgContactFilter.overCutoff.$Times.reciprocalMax.gfa.cluster`;

	$endhic_result_files .= "  $hicpro_matrix_file.CtgContactFilter.overCutoff.$Times.reciprocalMax.gfa.cluster";
}

print ENDHIC "perl $Bin/summarize_endhic_results.pl  --binsize $BinSize  --ctglencut $CtgSize_cutoff  --ctglenfile $contig_len_file   $endhic_result_files > z.EndHiC.$hicpro_matrix_file.$BinSize.results.summary.and.analysis  2> z.EndHiC.$hicpro_matrix_file.$BinSize.results.summary.and.analysis.cluster\n";
`perl $Bin/summarize_endhic_results.pl  --binsize $BinSize  --ctglencut $CtgSize_cutoff  --ctglenfile $contig_len_file   $endhic_result_files > z.EndHiC.$hicpro_matrix_file.$BinSize.results.summary.and.analysis  2> z.EndHiC.$hicpro_matrix_file.$BinSize.results.summary.and.analysis.cluster`;

close ENDHIC;


####################################################
################### Sub Routines ###################
####################################################
