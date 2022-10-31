#!/usr/bin/perl

=head1 Name

map_utg_to_ctg.pl  --  Mapping the hifiasm unitigs to contigs by shared reads and do statistics

=head1 Description

Input are Hifiasm gfa format files, output is a detailed table for mapping results

=head1 Version

  Author: Fan Wei, fanwei@caas.cn
  Version: 1.0,  Date: 2021/7/6
  Note:

=head1 Usage
  map_utg_to_ctg.pl <p_ctg_gfa_file>  <p_utg_gfa_file>
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

  perl ../map_utg_to_ctg.pl NB_hifiasm_l0.p_ctg.noseq.gfa NB_hifiasm_l0.p_utg.noseq.gfa > NB_hifiasm_l0.p_utg.noseq.gfa.mapped.p_ctg.list 2> NB_hifiasm_l0.p_utg.noseq.gfa.mapped.p_ctg.list.gfa
  
  perl ../map_utg_to_ctg.pl NB_hifiasm_l0.a_ctg.noseq.gfa NB_hifiasm_l0.p_utg.noseq.gfa > NB_hifiasm_l0.p_utg.noseq.gfa.mapped.a_ctg.list  2> NB_hifiasm_l0.p_utg.noseq.gfa.mapped.a_ctg.list.gfa

  Besides map p_utg to p_ctg and a_ctg, this program can also map r_utg to p_utg,  map r_utg to p_ctg and a_ctg.

  

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


my $ctg_gfa_file = shift;
my $utg_gfa_file = shift;

my %CtgReads;
my %UtgReads;
my %Utg;
my %UtgStrand;

my %CountCtg;
my %CountUtg;

my %CountCtgCov;
my %CountUtgCov;


my %CtgLen;
my %UtgLen;

my %UtgLinks;
my %CtgLinks;

my %Output;

my %UtgOrder;
my %UtgOrderStrand;


##parse the contig gfa file
open IN, $ctg_gfa_file;
while (<IN>) {
	
	#S       ptg000001l      *       LN:i:69740074   rd:i:73
	if(/^S\s+(\S+).+LN:i:(\d+)/){
		$CtgLen{$1} = $2;
		print STDERR $_;
	}

	##A       ptg000001l      41679   +       m64053_210226_012804/123404570/ccs      
	if(/^A\s+(\S+)\s+\d+\s+[+-]\s+(\S+)\s+\d+\s+(\d+)/){
		my $contig_id = $1;
		my $reads_id = $2;
		my $reads_len = $3;
		
		$CountCtg{$contig_id} ++;
		$CountCtgCov{$contig_id} += $reads_len;

		#print "$contig_id\t$reads_id\n";
		$CtgReads{$reads_id} = $contig_id;
	}

	if(/^L\s+(\S+)\s+[+-]\s+(\S+)/){
		my $ctg_id1 = $1;
		my $ctg_id2 = $2;
		$CtgLinks{$ctg_id1}{$ctg_id2} = 1;  ##double storing
		$CtgLinks{$ctg_id2}{$ctg_id1} = 1;  ##double storing
		print STDERR $_;
	}

}
close IN;


##parse the unitig gfa file
open IN, $utg_gfa_file;
while (<IN>) {
	
	if(/^S\s+(\S+).+LN:i:(\d+)/){
		$UtgLen{$1} = $2;
		print STDERR $_;
	}

	##A       ptg000001l      41679   +       m64053_210226_012804/123404570/ccs      
	if(/^A\s+(\S+)\s+\d+\s+([+-])\s+(\S+)\s+\d+\s+(\d+)/){
		my $unitig_id = $1;
		my $unitig_strand = $2;
		my $reads_id = $3;
		my $reads_len = $4;
		
		$CountUtg{$unitig_id} ++;
		$CountUtgCov{$unitig_id} += $reads_len;

		my $contig_id;
		if (exists $CtgReads{$reads_id}) {
			$contig_id = $CtgReads{$reads_id};
			$Utg{$unitig_id}{$contig_id} ++;
		}else{
			$contig_id = "None";
		}
		
		$UtgReads{$reads_id} = $unitig_id;
		$UtgStrand{$reads_id} = $unitig_strand;
	}
	
	#L       utg000004l      -       utg002489l      -       21379M  L1:i:23938496
	if(/^L\s+(\S+)\s+[+-]\s+(\S+)/){
		my $utg_id1 = $1;
		my $utg_id2 = $2;
		$UtgLinks{$utg_id1}{$utg_id2} = 1;  ##double storing
		$UtgLinks{$utg_id2}{$utg_id1} = 1;  ##double storing
		print STDERR $_;
	}

}
close IN;

#print Dumper \%Utg;



##Mapping the unitigs to contigs by shared reads_Ids, and store part of the output results
foreach my $utg_id (sort keys %Utg) {
	my $utg_p = $Utg{$utg_id};
	my @ids = keys %$utg_p;

	if(@ids == 1){
		my $ctg_id = $ids[0];
		my $shared_reads = $utg_p->{$ctg_id};
		my $ctg_reads = $CountCtg{$ctg_id};
		my $ctg_len = $CtgLen{$ctg_id};
		my $ctg_cov = sprintf("%.2f", $CountCtgCov{$ctg_id} / $ctg_len);

		my $utg_reads = $CountUtg{$utg_id};
		my $utg_len = $UtgLen{$utg_id};
		my $utg_cov = sprintf("%.2f", $CountUtgCov{$utg_id} / $utg_len);
		
		my $ratio_in_unitig = sprintf("%.4f", $shared_reads / $utg_reads);
		my $ratio_in_contig = sprintf("%.4f", $shared_reads / $ctg_reads);

		$Output{$ctg_id}{$utg_id} = "$utg_id\t$utg_len\t$utg_reads\t$utg_cov\t$ctg_id\t$ctg_len\t$ctg_reads\t$ctg_cov\t$utg_id\t$shared_reads\t$ratio_in_unitig\t$ratio_in_contig";

	}else{
		#print STDERR "Error:  unitig $utg_id mapped to more than one contigs\n";
	}
}


## pass the contig gfa file again, to sort the unitigs and get the mapping postion and strand information on the contigs
open IN, $ctg_gfa_file;
while (<IN>) {
	
	##A       ptg000001l      41679   +       m64053_210226_012804/123404570/ccs      
	if(/^A\s+(\S+)\s+(\d+)\s+([+-])\s+(\S+)\s+\d+\s+(\d+)/){
		my $contig_id = $1;
		my $contig_pos = $2;
		my $contig_strand = $3;
		my $reads_id = $4;
		if (!exists $UtgReads{$reads_id}) {
			#print STDERR "Error: $reads_id\t$contig_id\n";
			next;
		}
		my $unitig_id = $UtgReads{$reads_id};
		my $unitig_strand = $UtgStrand{$reads_id};
			
		if (! exists $UtgOrder{$contig_id}{$unitig_id}) {
			$UtgOrder{$contig_id}{$unitig_id} = $contig_pos;
			$UtgOrderStrand{$contig_id}{$unitig_id} = ($contig_strand eq $unitig_strand) ? "+": "-";
		}
		
	}

}
close IN;


##Final output of ordered unitigs for each contigs
print "UtgId\tUtgLength\tUtgReads#\tUtgDepth\tCtgId\tCtgLength\tCtgReads#\tCtgDepth\tUtgId\tSharedReads#\tRateInUtg\tRateInCtg\tBreak|Connect\tCtgId\t\tMapPosition\tMapStrand\tUtgNumInCtg\tUnique|repeat\n";

foreach my $contig_id (sort keys %Output) {
	
	my $order_p = $UtgOrder{$contig_id};
	
	my @ordered_unitigs = sort { $order_p->{$a} <=> $order_p->{$b} } keys %$order_p;
	
	##check the links between ordered neighboring unitigs, and find the broken site, which unitig do not have link with previous unitigs
	my $broken_status = 1; #whether the whole contigs has break point: 1 no; 0 yes;
	my %brokenIds;
	if (@ordered_unitigs >= 2) {
		for (my $i = 1; $i < @ordered_unitigs; $i ++) {
			my $current_utg = $ordered_unitigs[$i]; 
			
			##find the links between current and all previous unitigs
			my $link_status = 0;
			my $j = $i - 1;
			while($j >= 0){
				my $previous_utg = $ordered_unitigs[$j];
				if ( exists $UtgLinks{$current_utg}{$previous_utg}) {
					$link_status = 1;
					last;
				}
				$j --;
			}
			
			##If above failed, in addition, check the unmapped unitigs;
			my $previous_utg = $ordered_unitigs[$i - 1];
			if ($link_status == 0) {
				if (exists $UtgLinks{$current_utg}) {
					my $current_p = $UtgLinks{$current_utg};
					my @UnmappedUtgs = keys %$current_p;
					foreach my $unmapped_utg (@UnmappedUtgs) {
						if (exists $UtgLinks{$unmapped_utg}{$previous_utg}) {
							$link_status = 1;
							last;
						}

						##search another loop
						my $unmapped_p = $UtgLinks{$unmapped_utg};
						my @UnmappedUtgs2 = keys %$unmapped_p;
						foreach my $unmapped_utg2 (@UnmappedUtgs2) {
							if (exists $UtgLinks{$unmapped_utg2}{$previous_utg}) {
								$link_status = 1;
								last;
							}
						}


					}
				}
			}
			

			if ($link_status == 0) {
				$broken_status = 0;
				$brokenIds{$current_utg} = 1;
				#print STDERR "Error: $contig_id\t$current_utg\n";
			}
		}
	}

	my $unitig_count = @ordered_unitigs;
	
	##determine whether the contig is repeat by whether it has links to other contigs
	my $is_repeat = "Unique";
	if (exists $CtgLinks{$contig_id}) {
		$is_repeat = "Repeat";
		my $CtgLinks_p = $CtgLinks{$contig_id};
		my @CtgIds = keys %$CtgLinks_p;
		if (@CtgIds == 1 && $CtgIds[0] eq $contig_id) {
			$is_repeat = "Unique";
		}
	}

	##output each unitig information for each contigs
	foreach my $unitig_id (@ordered_unitigs) {
		
		my $contig_pos = $order_p->{$unitig_id};
		my $mapping_strand = $UtgOrderStrand{$contig_id}{$unitig_id};

		my $is_connected = (exists $brokenIds{$unitig_id} ) ? "Break" : "Connect";
		

		print $Output{$contig_id}{$unitig_id}."\t$is_connected\t$contig_id\t$contig_pos\t$mapping_strand\t$unitig_count\t$is_repeat\n";
	}

	##L       ptg000003l      +       ptg000186l      -       22731M  L1:i:86692131
	my $front_unitig = $ordered_unitigs[0];
	my $front_strand = ( $UtgOrderStrand{$contig_id}{$front_unitig} eq "+" ) ? "+" : "-";
	my $end_unitig = $ordered_unitigs[-1];
	my $end_strand = ( $UtgOrderStrand{$contig_id}{$end_unitig} eq "+" ) ? "-" : "+";
	

	print STDERR "L\t$contig_id\t+\t$front_unitig\t$front_strand\t0M\tL1:i:0\n";
	print STDERR "L\t$contig_id\t-\t$end_unitig\t$end_strand\t0M\tL1:i:0\n";
}



####################################################
################### Sub Routines ###################
####################################################
