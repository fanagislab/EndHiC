#!/usr/bin/perl

=head1 Name

endhic_ctgEnd_pipeline.pl  -- the pipeline of EndHic with a specified contig end size

=head1 Description

Invoke these programs sequentially: 
	(1) ctgContact_from_ctgEndContacts.pl   calculate the hic contact between all the contig ends with fixed length 
	(2) turningpoint_by_lineartransform.pl  get the turning point by adjusting the contacts data and linear transforming
	(3) scaffold_by_trueCtgContact.pl       scaffolding by link the contigs with true HiC contacts, output GFA file;
	(4) cluster_and_classify_GFA.pl         cluster from GFA and distinguish linear and circular topology
	(5) order_and_orient_GFA.pl             get order and orientation for each cluster from GFA
	(6) summerize_endhic_results.pl         summarize a set of Endhic results

=head1 Version

  Author: Fan Wei, fanwei@caas.cn
  Version: 1.0,  Date: 2021/7/23
  Note:

=head1 Usage
  
  endhic_ctgEnd_pipeline.pl [options] <contig_id_len_file>  <hicpro_bed_file> <hicpro_matrix_file> 
  --binsize <int>   bin size used in the hicpro result files, default=100000
  --binnum <int>    bin number used for contig end, contig end size = bin size x bin number, default=10
  --help            output help information to screen  

=head1 Example

  perl endhic_ctgEnd_pipeline.pl --binsize 100000 --binnum 10 contigs.id.len formal_100000_abs.bed formal_100000.matrix 
  perl endhic_ctgEnd_pipeline.pl --binsize 100000 --binnum 10 contigs.id.len formal_100000_abs.bed formal_1000000_iced.matrix

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;


##get options from command line into variables and set default values
my $Times;
my $BinSize;
my $BinNum;
my $expected_chromosome_number;
my $minimum_chromosome_length;
my $TimesStart;
my $TimesStop;
my $TimesStep;
my ($Verbose,$Help);
GetOptions(
	"binsize:i"=>\$BinSize,
	"binnum:i"=>\$BinNum,
	"chrnum:i"=>\$expected_chromosome_number,  ##hidden parameter, not for users
	"minchrlen:i"=>\$minimum_chromosome_length,  ##hidden parameter, not for users
	"mintimes:f"=>\$TimesStart,  ##hidden parameter, not for users
	"maxtimes:f"=>\$TimesStop,    ##hidden parameter, not for users
	"timesstep:f"=>\$TimesStep,    ##hidden parameter, not for users
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$BinSize ||= 100000;
$BinNum ||= 10;
$expected_chromosome_number ||= 0;
$minimum_chromosome_length ||= 10000000;
$TimesStart ||= 1.0;
$TimesStop  ||= 5.0;
$TimesStep  ||= 0.5;
die `pod2text $0` if (@ARGV == 0 || $Help);

my $contig_len_file = shift;
my $hicpro_bed_file = shift;
my $hicpro_matrix_file = shift;

my $endhic_result_files;

my $Contig_size_cutoff = $BinSize * $BinNum * 2;

print STDERR "Using Bin Size:  $BinSize;  Bin Number:  $BinNum\n";

my $raw_or_iced = ($hicpro_matrix_file =~ /_iced.matrix/) ? "iced" : "raw";

my $work_shell_file = "endhic.$BinSize.$BinNum.$raw_or_iced.sh";
open ENDHIC, ">$work_shell_file" || die "fail open $work_shell_file\n";

########################################################################################################


print ENDHIC "perl $Bin/ctgContact_from_ctgEndContacts.pl --binsize $BinSize --binnum $BinNum $contig_len_file $hicpro_bed_file $hicpro_matrix_file > $hicpro_matrix_file.$BinSize.$BinNum.CtgContact\n";
             `perl $Bin/ctgContact_from_ctgEndContacts.pl --binsize $BinSize --binnum $BinNum $contig_len_file $hicpro_bed_file $hicpro_matrix_file > $hicpro_matrix_file.$BinSize.$BinNum.CtgContact`;


if ($expected_chromosome_number > 0) {

	print ENDHIC "perl $Bin/ctgContact_normalize_distance.pl $hicpro_matrix_file.$BinSize.$BinNum.CtgContact > $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.distance\n";
				 `perl $Bin/ctgContact_normalize_distance.pl $hicpro_matrix_file.$BinSize.$BinNum.CtgContact > $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.distance`;

	print ENDHIC "perl $Bin/hcluster_contigs.pl --verbose -type min $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.distance $contig_len_file > $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.distance.hcluster.one 2> $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.distance.hcluster; rm $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.distance.hcluster.one\n";
				 `perl $Bin/hcluster_contigs.pl --verbose -type min $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.distance $contig_len_file > $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.distance.hcluster.one 2> $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.distance.hcluster; rm $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.distance.hcluster.one`;

	print ENDHIC "perl $Bin/hcluster_suitable_stop.pl --chr_num  $expected_chromosome_number --chr_len  $minimum_chromosome_length  $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.distance.hcluster  >  $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.distance.hcluster.need\n";
	             `perl $Bin/hcluster_suitable_stop.pl --chr_num  $expected_chromosome_number --chr_len  $minimum_chromosome_length  $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.distance.hcluster  >  $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.distance.hcluster.need`;
}

print ENDHIC "perl $Bin/turningpoint_by_lineartransform.pl $hicpro_matrix_file.$BinSize.$BinNum.CtgContact > $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.adjustTransform 2> $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.turningPoint\n";
             `perl $Bin/turningpoint_by_lineartransform.pl $hicpro_matrix_file.$BinSize.$BinNum.CtgContact > $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.adjustTransform 2> $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.turningPoint`;

my $raw_turning_point = `cat $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.turningPoint`;
$raw_turning_point = $1 if($raw_turning_point =~ /Turning point:\s+(\S+)/);
print STDERR "Contact turning point: $raw_turning_point\n";

my $raw_contact_cutoff;


for ($Times = $TimesStart; $Times <= $TimesStop; $Times += $TimesStep) {
	$Times = sprintf("%.1f",$Times);
	$raw_contact_cutoff = $Times * $raw_turning_point;
	print ENDHIC "perl $Bin/scaffold_by_trueCtgContact.pl   --contacts $raw_contact_cutoff $contig_len_file  $hicpro_matrix_file.$BinSize.$BinNum.CtgContact > $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.overCutoff.$Times.gfa \n";
	             `perl $Bin/scaffold_by_trueCtgContact.pl   --contacts $raw_contact_cutoff $contig_len_file  $hicpro_matrix_file.$BinSize.$BinNum.CtgContact > $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.overCutoff.$Times.gfa `;
	print ENDHIC "perl $Bin/scaffold_by_trueCtgContact.pl   --contacts $raw_contact_cutoff --reciprocalmax  $contig_len_file  $hicpro_matrix_file.$BinSize.$BinNum.CtgContact > $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.overCutoff.$Times.reciprocalMax.gfa \n";
	             `perl $Bin/scaffold_by_trueCtgContact.pl   --contacts $raw_contact_cutoff --reciprocalmax  $contig_len_file  $hicpro_matrix_file.$BinSize.$BinNum.CtgContact > $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.overCutoff.$Times.reciprocalMax.gfa `;
	
				  `perl $Bin/generate_GFAv1.2.pl $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.overCutoff.$Times.gfa > $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.overCutoff.$Times.gfa.v1.2.gfa`;
				  `perl $Bin/generate_GFAv1.2.pl $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.overCutoff.$Times.reciprocalMax.gfa > $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.overCutoff.$Times.reciprocalMax.gfa.v1.2.gfa`;
	
	print ENDHIC "perl $Bin/cluster_and_classify_GFA.pl  $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.overCutoff.$Times.reciprocalMax.gfa > $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.overCutoff.$Times.reciprocalMax.gfa.topology\n";
	             `perl $Bin/cluster_and_classify_GFA.pl  $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.overCutoff.$Times.reciprocalMax.gfa > $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.overCutoff.$Times.reciprocalMax.gfa.topology`;
	
	print ENDHIC "perl $Bin/order_and_orient_GFA.pl --size $Contig_size_cutoff $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.overCutoff.$Times.reciprocalMax.gfa $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.overCutoff.$Times.reciprocalMax.gfa.topology > $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.overCutoff.$Times.reciprocalMax.gfa.cluster\n";
	             `perl $Bin/order_and_orient_GFA.pl --size $Contig_size_cutoff $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.overCutoff.$Times.reciprocalMax.gfa $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.overCutoff.$Times.reciprocalMax.gfa.topology > $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.overCutoff.$Times.reciprocalMax.gfa.cluster`;
	
	$endhic_result_files .= "  $hicpro_matrix_file.$BinSize.$BinNum.CtgContact.overCutoff.$Times.reciprocalMax.gfa.cluster";
}


print ENDHIC "perl $Bin/summarize_endhic_results.pl  --binsize $BinSize --ctglencut $Contig_size_cutoff  --ctglenfile $contig_len_file  $endhic_result_files  > $hicpro_matrix_file.$BinSize.$BinNum.results.summary  2> $hicpro_matrix_file.$BinSize.$BinNum.results.summary.cluster\n";
             `perl $Bin/summarize_endhic_results.pl  --binsize $BinSize --ctglencut $Contig_size_cutoff  --ctglenfile $contig_len_file  $endhic_result_files  > $hicpro_matrix_file.$BinSize.$BinNum.results.summary  2> $hicpro_matrix_file.$BinSize.$BinNum.results.summary.cluster`;

close ENDHIC;


####################################################
################### Sub Routines ###################
####################################################
