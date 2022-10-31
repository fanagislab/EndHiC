#!/usr/bin/perl

=head1 Name

revise_iced_matrix.pl   -- revise the hic-pro iced_matrix file

=head1 Description

In the hic-pro iced_matrix file, the normalization of some contacts may cause huge mistake;
For example, the raw contact value is very low, but the normalized contact value is quite high
So we revise the hic-pro iced_matrix by replacing the normalized contact with raw contact, whose
raw contact value is less than a given cutoff.

=head1 Version

  Author: Fan Wei, fanw@caas.cn
  Version: 1.0,  Date: 2021/8/18
  Note:

=head1 Usage
  
  revise_iced_matrix.pl [options] <raw_matrix_file> <iced_matrix_file>
  --rawcontacts <int>  cutoff for the raw contact, default=10 
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

  revise_iced_matrix.pl --rawcontacts 10 formal_100000.matrix formal_100000_iced.matrix > formal_100000_iced.matrix.revised

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
$recover_cutoff ||= 10;
die `pod2text $0` if (@ARGV == 0 || $Help);

my $raw_matrix_file = shift;
my $ice_matrix_file = shift;

my %RawMatrix;



open IN, $raw_matrix_file || die "fail open $raw_matrix_file";
while (<IN>) {
	chomp;
	my ($bin1, $bin2, $contacts) = split /\s+/;
	$RawMatrix{$bin1}{$bin2} = $contacts;
}
close IN;

open IN, $ice_matrix_file  || die "fail open $ice_matrix_file ";
while (<IN>) {
	chomp;
	my ($bin1, $bin2, $IcedContacts) = split /\s+/;
	my $rawContacts = $RawMatrix{$bin1}{$bin2};
	##print "$bin1\t$bin2\t$rawContacts\t$IcedContacts\n";
	$IcedContacts = $rawContacts if($rawContacts <= $recover_cutoff);
	$IcedContacts = int($IcedContacts);
	print "$bin1\t$bin2\t$IcedContacts\t$rawContacts\n";
}
close IN;


####################################################
################### Sub Routines ###################
####################################################
