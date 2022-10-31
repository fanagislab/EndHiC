#!/usr/bin/perl

=head1 Name

ctgContact_from_ctgEndContacts.pl  -- calculate HiC contact for both contig ends: the head and the tail

=head1 Description

This program calculated the hic contact counts between all the contig heads and tails, whose length are 
defined by a specified number of 100-kb bins. The recommonded (default) bin number is 10, then both the 
contig head and contig tail length used are 1 Mb.

=head1 Version

  Author: Fan Wei, fanwei@caas.cn
  Version: 1.0,  Date: 2021/7/6
  Note:

=head1 Usage
  ctgContact_from_ctgEndContacts.pl <contigs_used.len>  <formal_100000_abs.bed>  <formal_100000.matrix>
  --binsize <int>  size of the input Bins, unit in bp, default=100000
  --binnum <int>  number of bins included in contig head or contig tail, -1 for half-contig bins, default=10
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

  perl ../ctgContact_from_ctgEndContacts.pl --binsize 100000 --binnum 10 contigs_used.len formal_100000_abs.bed formal_100000.matrix > formal_100000.matrix.CtgContact
  

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;


##get options from command line into variables and set default values
my $HicPro_BinSize;
my $HicPro_BinNum;
my ($Verbose,$Help);
GetOptions(
	"binsize:i"=>\$HicPro_BinSize,
	"binnum:i"=>\$HicPro_BinNum,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$HicPro_BinSize ||= 100000;
$HicPro_BinNum ||= 10;
die `pod2text $0` if (@ARGV == 0 || $Help);

my $contig_end_size = ($HicPro_BinNum > 0) ? $HicPro_BinNum * $HicPro_BinSize :  $HicPro_BinSize;

if ($HicPro_BinNum > 0) {
	print STDERR "Size of contig head or tail end: $HicPro_BinNum x $HicPro_BinSize = $contig_end_size bp\n";
}else{
	print STDERR "Size of contig head or tail end: half of the contig size, minimum is $contig_end_size\n";
}

my $Contig_used_file = shift;
my $Contig_bin_file = shift;
my $Bin_matrix_file = shift;

my %UsedContigs;
my %Bin2Ctg;

my %ContactData;

##only do statistics for contigs included in $Contig_used_file
open IN, $Contig_used_file || die "fail open $Contig_used_file\n";
while (<IN>) {
	
	if(/(\S+)\s+(\d+)/){
		if($2 < $contig_end_size*2){
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
my %BinData;
while (<IN>) {
	#ptg000001l      0       1000000 1
	my ($CtgId, $BinId) = ($1,$2) if(/^(\S+)\s+\d+\s+\d+\s+(\d+)/);
	next if(! exists $UsedContigs{$CtgId});
	
	push @{$BinData{$CtgId}}, $BinId;
}
close IN;
foreach my $ctgId (sort keys %BinData) {
	my $ctg_p = $BinData{$ctgId};
	my $ctg_len = $UsedContigs{$ctgId};
	

	my $Used_BinNum = ($HicPro_BinNum > 0) ? $HicPro_BinNum : int($ctg_len/$HicPro_BinSize/2) ;
	
	for (my $i = 0; $i<$Used_BinNum; $i ++) {
		my $BinId = $ctg_p->[$i];
		$Bin2Ctg{$BinId} = $ctgId."-head";
		#print $BinId,"\t",$ctgId."-head\n";

	}
	for (my $i = 1; $i<$Used_BinNum+1; $i ++) {
		my $BinId = $ctg_p->[-$i];
		$Bin2Ctg{$BinId} = $ctgId."-tail" if(!exists $Bin2Ctg{$BinId});
		#print $BinId,"\t",$ctgId."-tail\n";
	}
}

#print STDERR Dumper \%Bin2Ctg;


##this is hic-pro generated contact map matrix file
open IN, $Bin_matrix_file || die "fail open $Bin_matrix_file\n";
while (<IN>) {
	#1       1       9.183982
	my ($BinId1, $BinId2, $Contact) = ($1, $2, $3) if(/^(\d+)\s+(\d+)\s+(\S+)/);
	next if(! exists $Bin2Ctg{$BinId1}  || ! exists $Bin2Ctg{$BinId2});

	my $CtgId1 = $Bin2Ctg{$BinId1};
	my $CtgId2 = $Bin2Ctg{$BinId2};
	
	#print STDERR "$CtgId1\t$CtgId2\n";

	my $CtgId1_str = $CtgId1;
	$CtgId1_str =~ s/-head//;
	$CtgId1_str =~ s/-tail//;
	my $CtgId2_str = $CtgId2;
	$CtgId2_str =~ s/-head//;
	$CtgId2_str =~ s/-tail//;

	$ContactData{$CtgId1}{$CtgId2} += $Contact if($CtgId1_str ne $CtgId2_str);

}
close IN;

#print Dumper \%ContactData;


##output the contact density between contig pairs
my @Output;
foreach  my $CtgId1 (sort keys %ContactData) {
	my $CtgId1_p = $ContactData{$CtgId1};
	foreach my $CtgId2 (sort keys %$CtgId1_p) {
		
		my ($CtgId1_str ,$CtgId1_pos) = ($1,$2) if($CtgId1 =~ /(\w+)-(\w+)/);
		my ($CtgId2_str ,$CtgId2_pos) = ($1,$2) if($CtgId2 =~ /(\w+)-(\w+)/);
		
		my $Used_BinNum1 = ($HicPro_BinNum > 0) ? $HicPro_BinNum : int($UsedContigs{$CtgId1_str}/$HicPro_BinSize/2) ;
		my $Used_BinNum2 = ($HicPro_BinNum > 0) ? $HicPro_BinNum : int($UsedContigs{$CtgId2_str}/$HicPro_BinSize/2) ;
		
		my $total = $CtgId1_p->{$CtgId2};
		#$total = $total / ($Used_BinNum1 * $Used_BinNum2);  #normalized by bin number
		
		$total = sprintf("%.2f", $total);

		#print "$CtgId1_str\t$CtgId2_str\t$total\t$CtgId1_pos\t$CtgId2_pos\n"; ##output the total contacts between two DNA fragments
		push @Output, [$CtgId1_str,$CtgId2_str,$total,$CtgId1_pos,$CtgId2_pos, $Used_BinNum1, $Used_BinNum2];
	}

}

@Output = sort {$b->[2] <=> $a->[2]} @Output;

print "#CtgId1\tCtgId2\tEndContact\tCtg1Pos\tCtg2Pos\tUsedBinNum1\tUsedBinNum2\n";
foreach my $p (@Output) {
	my $line = join("\t",@$p);
	print $line."\n";
}

