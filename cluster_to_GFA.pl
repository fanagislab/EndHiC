#!/usr/bin/perl

=head1 Name

cluster_to_GFA.pl  -- convert cluster format to GFA format

=head1 Description


=head1 Version

  Author: Fan Wei, fanwei@caas.cn
  Version: 1.0,  Date: 2021/7/23
  Note:

=head1 Usage
  
  cluster_to_GFA.pl  <ctg_id_len_file>  <cluster_file>
  --help            output help information to screen  

=head1 Example

cluster_to_GFA.pl contigs.id.len z.EndHiC.A.results.summary.cluster > z.EndHiC.A.results.summary.cluster.gfa

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


my $Contig_used_file = shift;
my $cluster_file = shift;

my %CtgLen;

my %Ctgs;


##get contig id and length data
open IN, $Contig_used_file || die "fail open $Contig_used_file\n";
while (<IN>) {
	if(/(\S+)\s+(\d+)/){
		$CtgLen{$1} = $2;
	}
}
close IN;


##S       ptg000007l      *       LN:i:180134246
##L       ptg000003l      +       ptg000004l      -       0M      L1:i:0  RM:i:1  Contact:i:2603.00


my $S_out;
my $L_out;

open IN, $cluster_file || die "fail open $cluster_file";
while (<IN>) {
	next if(/^\#/ || /Unclustered/);
	chomp;
	my ($cluster_id, $ctg_num, $cluster_len, $robustness, $ctg_str) = split /\s+/;
	
	my @t = split /;/, $ctg_str;
	
	for (my $i=0; $i<@t; $i++) {
		my ($ctgid, $ctgstrand) = ($1,$2) if($t[$i] =~ /(\w+)([+-])/);
		my $ctglen = $CtgLen{$ctgid};
		$S_out .= "S\t$ctgid\t*\tLN:i:$ctglen\n";
	}
	
	for (my $i=0; $i<@t - 1; $i++) {
		my ($ctg1id, $ctg1strand) = ($1,$2) if($t[$i] =~ /(\w+)([+-])/);
		my ($ctg2id, $ctg2strand) = ($1,$2) if($t[$i+1] =~ /(\w+)([+-])/);
		$L_out .= "L\t$ctg1id\t$ctg1strand\t$ctg2id\t$ctg2strand\t0M\n";
	}
		
	
}
close IN;


print $S_out;
print $L_out;

