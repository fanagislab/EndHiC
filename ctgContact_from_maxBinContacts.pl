#!/usr/bin/perl

=head1 Name

ctgContact_from_maxBinContacts.pl  -- calculate contig contacts from the max bin contacts 

=head1 Description

This program calculated the hic contact counts between all the contigs, by using the max contact count
of the bins belong to each compared contigs.

Take the hic-pro bin contact matrix as input data, and output a tabular format file

To balance the effectiveness and robustness, bin size 1000-Kb are recommonded.

Contigs with length shorter than 2 times of the bin size are filtered out.

=head1 Version

  Author: Fan Wei, fanwei@caas.cn
  Version: 1.0,  Date: 2021/7/6
  Note:

=head1 Usage
  ctgContact_from_maxBinContacts.pl <contig_lenght_file>  <hic_pro_bed_file> <hic_pro_matrix_file>
  --binsize <int>  size of the input Bins, unit in bp, default=1000000
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

  perl ../ctgContact_from_maxBinContacts.pl --binsize 1000000 contigs_used.len formal_1000000_abs.bed  formal_1000000.matrix > formal_1000000.matrix.ctgcontact
  
  perl ../ctgContact_from_maxBinContacts.pl --binsize 1000000 contigs_used.len formal_1000000_abs.bed  formal_1000000_iced.matrix > formal_1000000_iced.matrix.ctgcontact


=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;


##get options from command line into variables and set default values
my $HicPro_BinSize;
my ($Verbose,$Help);
GetOptions(
	"binsize:i"=>\$HicPro_BinSize,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$HicPro_BinSize ||= 1000000;
die `pod2text $0` if (@ARGV == 0 || $Help);

my $Contig_used_file = shift;
my $Contig_bin_file = shift;
my $Bin_matrix_file = shift;

my %UsedContigs;
my %Bin2Ctg;
my %CtgBinNum;
my %CtgBinStart;

my %ContactData;
my %Bin1Data;
my %Bin2Data;



##only do statistics for contigs included in $Contig_used_file
open IN, $Contig_used_file || die "fail open $Contig_used_file\n";
while (<IN>) {
	if(/(\S+)\s+(\d+)/){
		if($2 < $HicPro_BinSize*2){
			#print STDERR "$1\tis shorter the two times of contig end size, and be ignored\n";
			next;
		}
		$UsedContigs{$1} = $2;
	}
}
close IN;

##print Dumper \%UsedContigs;

##this is hic-pro generated .bed file
open IN, $Contig_bin_file || die "fail open $Contig_bin_file\n";
while (<IN>) {
	#ptg000001l      0       1000000 1
	my ($CtgId, $BinId) = ($1,$2) if(/^(\S+)\s+\d+\s+\d+\s+(\d+)/);
	next if(! exists $UsedContigs{$CtgId});

	$Bin2Ctg{$BinId} = $CtgId ;
	$CtgBinNum{$CtgId} = $BinId;
	if (! exists $CtgBinStart{$CtgId}) {
		$CtgBinStart{$CtgId} = $BinId;
	}
}
close IN;

##print Dumper \%Bin2Ctg;

#print Dumper \%CtgBinNum;

##this is hic-pro generated contact map matrix file
open IN, $Bin_matrix_file || die "fail open $Bin_matrix_file\n";
while (<IN>) {
	#1       1       9.183982
	my ($BinId1, $BinId2, $Contact) = ($1, $2, $3) if(/^(\d+)\s+(\d+)\s+(\S+)/);
	next if(! exists $Bin2Ctg{$BinId1}  || ! exists $Bin2Ctg{$BinId2});

	my $CtgId1 = $Bin2Ctg{$BinId1};
	my $CtgId2 = $Bin2Ctg{$BinId2};

	if($CtgId1 ne $CtgId2){
		push @{$ContactData{$CtgId1}{$CtgId2}}, $Contact;
		push @{$Bin1Data{$CtgId1}{$CtgId2}}, $BinId1;
		push @{$Bin2Data{$CtgId1}{$CtgId2}}, $BinId2;
	}

}
close IN;

#print Dumper \%ContactData;

##output the contact density between contig pairs
my @Output;
foreach  my $CtgId1 (sort keys %ContactData) {
	my $CtgId1_p = $ContactData{$CtgId1};
	foreach my $CtgId2 (sort keys %$CtgId1_p) {
		my $contact_p = $CtgId1_p->{$CtgId2};
		my $Bin1_p = $Bin1Data{$CtgId1}{$CtgId2};
		my $Bin2_p = $Bin2Data{$CtgId1}{$CtgId2};
		
		my $ctg1_binNum = $CtgBinNum{$CtgId1};
		my $ctg1_binStart = $CtgBinStart{$CtgId1};
		my $ctg2_binNum = $CtgBinNum{$CtgId2};
		my $ctg2_binStart = $CtgBinStart{$CtgId2};


		my $num = @$contact_p;
		my $total = 0;
		my $avg = 0;
		my $max = 0;
		my $max_Bin1 = $ctg1_binStart;
		my $max_Bin2 = $ctg2_binStart;
		my $i = 0;
		foreach my $this_contact (@$contact_p) {
			if ($this_contact > $max) {
				$max = $this_contact;
				$max_Bin1 = $Bin1_p->[$i];
				$max_Bin2 = $Bin2_p->[$i];
			}
			$total += $this_contact;
			$i++;
		}
		$avg = $total / $num;
		
		$total = sprintf("%.2f", $total);
		$avg = sprintf("%.2f", $avg);
		$max = sprintf("%.2f", $max);

		my $distance1;
		my $postion_rate1;
		my $headtail_status1 = ($max_Bin1 - $ctg1_binStart <= $ctg1_binNum - $max_Bin1) ? "Head" : "Tail";
		if ($headtail_status1 eq "Head") {
			$postion_rate1 = ($max_Bin1 - $ctg1_binStart) / ($ctg1_binNum -  $ctg1_binStart + 1);
			$distance1 =  $max_Bin1 - $ctg1_binStart;
		}else{
			$postion_rate1 = ($max_Bin1 - $ctg1_binStart + 1) / ($ctg1_binNum -  $ctg1_binStart + 1);
			$distance1 = $ctg1_binNum - $max_Bin1;
		}
		$postion_rate1 = sprintf("%.2f",$postion_rate1);
		
		my $distance2;
		my $postion_rate2;
		my $headtail_status2 = ($max_Bin2 - $ctg2_binStart <= $ctg2_binNum - $max_Bin2) ? "Head" : "Tail";
		if ($headtail_status2 eq "Head") {
			$postion_rate2 = ($max_Bin2 - $ctg2_binStart) / ($ctg2_binNum -  $ctg2_binStart + 1);
			$distance2 = $max_Bin2 - $ctg2_binStart;
		}else{
			$postion_rate2 = ($max_Bin2 - $ctg2_binStart + 1) / ($ctg2_binNum -  $ctg2_binStart + 1);
			$distance2 = $ctg2_binNum - $max_Bin2;
		}
		$postion_rate2 = sprintf("%.2f",$postion_rate2);
		
		#print "$CtgId1\t$CtgId2\t$max\t$headtail_status1\t$headtail_status2\t$postion_rate1\t$postion_rate2\t$distance1\t$distance2\t";
		#print "$CtgId1\t$ctg1_binStart-$ctg1_binNum\t$max_Bin1\t";
		#print "$CtgId2\t$ctg2_binStart-$ctg2_binNum\t$max_Bin2\n";
		push @Output, [$CtgId1,$CtgId2,$max,$headtail_status1,$headtail_status2,$postion_rate1,$postion_rate2,$distance1,$distance2,$CtgId1,"$ctg1_binStart-$ctg1_binNum",$max_Bin1,$CtgId2,"$ctg2_binStart-$ctg2_binNum",$max_Bin2, $avg];
	
	}

}

@Output = sort {$b->[2] <=> $a->[2]} @Output;

print "#CtgId1\tCtgId2\tMaxContact\tCtg1Pos\tCtg2Pos\tCtg1PosRate\tCtg2PosRate\tDist1ToEnd\tDist2ToEnd\tCtgId1\tCtg1BinRange\tCtg1Bin\tCtgId2\tCtg2BinRange\tCtg2Bin\tAvgContact\n";
foreach my $p (@Output) {
	my $line = join("\t",@$p);
	print $line."\n";
}
