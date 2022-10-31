#!/usr/bin/perl

=head1 Name

cluster_compare.pl  -- compare two cluster format files

=head1 Description 


=head1 Version

  Author: Fan Wei, fanwei@caas.cn
  Version: 1.0,  Date: 2021/7/23
  Note:

=head1 Usage
  
  cluster_compare.pl  <cluster_file1> <cluster_file2>

  --help            output help information to screen  

=head1 Example


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

my $agp_file = shift;

my %Info;

open IN, $agp_file || die "fail";
while (<IN>) {
	chomp;
	my @t = split /\t/;
	next if($t[4] ne "W");

	my $chr_id = $t[0];
	my $ctg_id = $t[5];
	my $ctg_len = $t[7];
	my $strand = $t[8];

	$Info{$chr_id}{"count"} ++;
	$Info{$chr_id}{"length"} += $ctg_len;
	$Info{$chr_id}{"ctgs_str"} .= "$ctg_id$strand;";
}
close IN;

print "#Cluster_id	contig_count	cluster_length	robustness[max:90]	contigs_order_orientation\n";
foreach my $chr_id (sort keys %Info) {
	my $chr_p = $Info{$chr_id};
	my $count = $chr_p->{"count"};
	my $length = $chr_p->{"length"};
	my $ctgs_str = $chr_p->{"ctgs_str"};
	$ctgs_str =~ s/;$//;
	print "$chr_id\t$count\t$length\t*\t$ctgs_str\n";
}