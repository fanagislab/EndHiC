#!/usr/bin/perl

=head1 Name

order_and_orient_GFA.pl  --  cluster, order and orient contigs from the GFA format file

=head1 Description

Take the result of scaffold_by_trueCtgContact.pl and cluster_and_classify_GFA.pl as input data

The output is in format: Cluster_id     ctg_count       cluster_len     ctgs_order_orientation

Before the major work, the circular cluster is broken by invoking "linearize_GFA.pl"

=head1 Version

  Author: Fan Wei, fanw@caas.cn
  Version: 1.0,  Date: 2021/8/22
  Note:

=head1 Usage
  
  order_and_orient_GFA.pl <gfa_file> <gfa_cluster_file>
  --size <int>  size cutoff of the cluster, default=0
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

   order_and_orient_GFA.pl --size 6000000 formal_100000_iced.matrix.revised.100000.30.CtgContact.overCutoff.5.0.reciprocalMax.gfa formal_100000_iced.matrix.revised.100000.30.CtgContact.overCutoff.5.0.reciprocalMax.gfa.cluster > formal_100000_iced.matrix.revised.100000.30.CtgContact.overCutoff.5.0.reciprocalMax.gfa.cluster.order.orient

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;


##get options from command line into variables and set default values
my $Size_cutoff;
my ($Verbose,$Help);
GetOptions(
	"size:i"=>\$Size_cutoff,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Size_cutoff ||= 0;
die `pod2text $0` if (@ARGV == 0 || $Help);

my $gfa_file = shift;
my $gfa_cluster_file = shift;

my $linearized_gfa_file = "$gfa_file.linear";
`perl $Bin/linearize_GFA.pl $gfa_file $gfa_cluster_file > $linearized_gfa_file`;



my %CtgLen;
my %Links;
my @EndNodes;
my %OrderData;

my @Cluster;

##this is the linearized gfa file in which no branches existing
open IN, $linearized_gfa_file || die "fail open $linearized_gfa_file";
while (<IN>) {
	
	my $line = $_;

	if (/^S/) {
		##S       ptg000008l      *       LN:i:61985296
		my ($ctgId, $ctgLen) = ($1,$2) if(/^S\s+(\S+).+LN:i:(\d+)/);
		$CtgLen{$ctgId} = $ctgLen;
	}
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

		push @{$Links{$ctg1}}, $ctg2;
		push @{$Links{$ctg2}}, $ctg1;
	}
}
close IN;

#print Dumper \%OrderData;;

##Find the end contigs which has one and only one link to other contigs 
foreach my $ctg (sort keys %Links) {
	my $ctp_p = $Links{$ctg};
	if(@$ctp_p == 1){
		#print $ctg."\n";
		push @EndNodes, $ctg;
	}
}

#print Dumper \%EndNodes;


##prepare the order data and make information in the reverse and complementary
open IN, $linearized_gfa_file || die "fail open $linearized_gfa_file";
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
		
		$OrderData{$ctg1.$ctg1_strand} = $ctg2.$ctg2_strand;
		
		my $new_line;
		if ($ctg1_strand ne $ctg2_strand) {
			#$new_line = "$type\t$ctg2\t$ctg1_strand\t$ctg1\t$ctg2_strand\t$match\t$L1\t$RM\n";
			$OrderData{$ctg2.$ctg1_strand} = $ctg1.$ctg2_strand;
		}else{
			my $new_strand = ($ctg1_strand eq "+") ? "-" : "+";
			#$new_line = "$type\t$ctg2\t$new_strand\t$ctg1\t$new_strand\t$match\t$L1\t$RM\n";
			$OrderData{$ctg2.$new_strand} = $ctg1.$new_strand;
		}
	}
}
close IN;

##print Dumper \%OrderData;

##generate clusters without no redundance
my %PassedEndNodes;
foreach my $start_ctg (@EndNodes) {
	
	next if (exists $PassedEndNodes{$start_ctg});

	my $this_ctg = (exists $OrderData{$start_ctg."+"}) ? $start_ctg."+" : $start_ctg."-";
	my @ary;
	push @ary, $this_ctg;
	while (1) {
		if (exists $OrderData{$this_ctg}) {
			my $next_ctg = $OrderData{$this_ctg};
			push @ary, $next_ctg;
			$this_ctg = $next_ctg;
			
		}else{
			last;
		}
	}
	my $stop_ctg = $1 if($ary[-1] =~ /^(\w+)[+-]/);
	push @Cluster, \@ary;
	
	$PassedEndNodes{$start_ctg} = 1;
	$PassedEndNodes{$stop_ctg} = 1;

}

#print Dumper \@Cluster;


##get linked contigs
my %LinkedCtgs;
foreach my $p (@Cluster) {
	foreach my $str (@$p) {
		my $ctgId = $1 if($str =~ /^(\w+)[+-]/);
		$LinkedCtgs{$ctgId} = 1;
	}
}

#print Dumper \%LinkedCtgs;

#Add un-linked contigs into @Cluster
foreach my $ctgId (keys %CtgLen) {
	if (! exists $LinkedCtgs{$ctgId}) {
		push @Cluster, [$ctgId."+"];
	}
}

##print Dumper \@Cluster;

##calculate the total length of each cluster
my %ClusterLen;
for (my $i=0; $i<@Cluster; $i++) {
	my $p = $Cluster[$i];
	my $cluster_len = 0;
	foreach my $str (@$p) {
		my $ctgId = $1 if($str =~ /^(\w+)[+-]/);
		my $ctgLen = $CtgLen{$ctgId};
		$cluster_len += $ctgLen;
	}
	$ClusterLen{$i} = $cluster_len;
}

#print Dumper \%ClusterLen;


##output the final result
print "#Cluster_id\tctg_count\tcluster_len\trobustness\tctgs_order_orientation\n";
my $loop = "01";
foreach my $i (sort {$ClusterLen{$b} <=> $ClusterLen{$a}} keys %ClusterLen) {
	
	my $cluster_id = "Cluster_$loop";
	my $cluster_len = $ClusterLen{$i};
	my $cluster_p = $Cluster[$i];
	my $ctg_count = @$cluster_p;
	my $ctg_str = join(";", @$cluster_p);
	print "$cluster_id\t$ctg_count\t$cluster_len\t1\t$ctg_str";
	print "\t[Type:Unclustered]" if($ctg_count == 1 && $cluster_len < $Size_cutoff);
	print "\n";

	$loop ++;
}


`rm $linearized_gfa_file`;

####################################################
################### Sub Routines ###################
####################################################
