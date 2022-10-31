#!/usr/bin/perl

=head1 Name



=head1 Description



=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage

  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple



=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;


##get options from command line into variables and set default values
my $recover_cutoff;
my ($Verbose,$Help);
GetOptions(
	"rawcontacts:i"=>\$recover_cutoff,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$recover_cutoff ||= 100;
die `pod2text $0` if (@ARGV == 0 || $Help);

my $Contig_bin_file = shift;
my $raw_matrix_file = shift;
my $ice_matrix_file = shift;

my %Bin2Ctg;
my %IcedMatrix;


##this is hic-pro generated .bed file
open IN, $Contig_bin_file || die "fail open $Contig_bin_file\n";
while (<IN>) {
	#ptg000001l      0       1000000 1
	my ($CtgId, $BinId) = ($1,$2) if(/^(\S+)\s+\d+\s+\d+\s+(\d+)/);

	$Bin2Ctg{$BinId} = $CtgId ;
	
}
close IN;



open IN, $ice_matrix_file || die "fail open $ice_matrix_file";
while (<IN>) {
	chomp;
	my ($bin1, $bin2, $contacts) = split /\s+/;
	$IcedMatrix{$bin1}{$bin2} = $contacts;
}
close IN;


print "ctg1\tctg2\tbin1\tbin2\trawContacts\ticedContacts\n";
open IN, $raw_matrix_file  || die "fail open $raw_matrix_file";
while (<IN>) {
	chomp;
	my ($bin1, $bin2, $rawContacts) = split /\s+/;

	my $icedContacts = (exists $IcedMatrix{$bin1}{$bin2}) ? $IcedMatrix{$bin1}{$bin2} : 0;
	
	my $ctg1 = $Bin2Ctg{$bin1};
	my $ctg2 = $Bin2Ctg{$bin2};

	print "$ctg1\t$ctg2\t$bin1\t$bin2\t$rawContacts\t$icedContacts\n";
}
close IN;


####################################################
################### Sub Routines ###################
####################################################
