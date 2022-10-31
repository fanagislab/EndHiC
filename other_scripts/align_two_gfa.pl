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
my $gfa1_mark;
my $gfa2_mark;
my ($Verbose,$Help);
GetOptions(
	"mark:s"=>\$gfa1_mark,
	"mark:s"=>\$gfa2_mark,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$gfa1_mark ||= "A";
$gfa2_mark ||= "B";
die `pod2text $0` if (@ARGV == 0 || $Help);


my $gfa1_file = shift;
my $gfa2_file = shift;

my @Contigs;

##this is the linearized gfa file in which no branches existing
open IN, $gfa1_file || die "fail open $gfa1_file";
while (<IN>) {
	
	chomp;
	my @t = split /\s+/;

	if (/^S/) {
		push @Contigs, $t[1];
		$t[1] = $t[1].".".$gfa1_mark;

	}
	if (/^L/) {
		#L       ptg000002l      +       ptg000026l      +       0M      L1:i:0  RM:i:1
		$t[1] = $t[1].".".$gfa1_mark;
		$t[3] = $t[3].".".$gfa1_mark;
	}

	my $line = join("\t", @t);
	print $line."\n";
}
close IN;




##this is the linearized gfa file in which no branches existing
open IN, $gfa2_file || die "fail open $gfa2_file";
while (<IN>) {
	
	chomp;
	my @t = split /\s+/;

	if (/^S/) {
		$t[1] = $t[1].".".$gfa2_mark;

	}
	if (/^L/) {
		#L       ptg000002l      +       ptg000026l      +       0M      L1:i:0  RM:i:1
		$t[1] = $t[1].".".$gfa2_mark;
		$t[3] = $t[3].".".$gfa2_mark;
	}

	my $line = join("\t", @t);
	print $line."\n";
}
close IN;


##L       ptg000027l      +       ptg000037l      -       0M      L1:i:0  RM:i:1
foreach my $ctgId (@Contigs) {
	my $ctgIdMark1 = $ctgId.".".$gfa1_mark;
	my $ctgIdMark2 = $ctgId.".".$gfa2_mark;
	print "L\t$ctgIdMark1\t+\t$ctgIdMark2\t+\t0M\tAdd:1\n";
	print "L\t$ctgIdMark1\t-\t$ctgIdMark2\t-\t0M\tAdd:1\n";
}

####################################################
################### Sub Routines ###################
####################################################
