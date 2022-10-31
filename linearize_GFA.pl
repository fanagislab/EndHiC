#!/usr/bin/perl

=head1 Name

linear_GFA.pl  --  linearize the GFA file

=head1 Description

Break the most lowest link in a circular cluster

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage
  
  linearize_GFA.pl  <gfa_file>  <gfa_cluster_file>
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

   linearize_GFA.pl gfa_file gfa_cluster_file > linearized_gfa_file

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

my $gfa_file = shift;
my $gfa_cluster_file = shift;

my %LinkData;
my %DeletedLines;

open IN, $gfa_file || die "fail open $gfa_file";
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
		my $contact = $1 if($t[8] =~ /Contact:i:(\d+)/);

		$LinkData{$line} = $contact;
	}
}
close IN;

open IN, $gfa_cluster_file || die "fail open $gfa_cluster_file";
while (<IN>) {
	my @t = split /\s+/;
	
	my $cluster_id = shift @t;
	my $ctgs_num = shift @t;
	my $type = pop @t;
	my @ctgs = @t;

	

	if ($type =~ /Circular/) {
		#print $cluster_id."\n";
		my %LinesInCluster;
		foreach my $ctg (@ctgs) {
			foreach my $line (keys %LinkData) {
				my $contact = $LinkData{$line};
				if ($line =~ /$ctg/) {
					$LinesInCluster{$line} = $contact;
				}
			}
		}
		
		my $min_contact = 1000000000;
		my $min_line;
		foreach my $line (keys %LinesInCluster) {
			my $contact = $LinesInCluster{$line};
			if ($contact < $min_contact) {
				$min_contact = $contact;
				$min_line = $line;
			}
		}
		#print $min_line."\n";
		$DeletedLines{$min_line} = 1;

	}

}	
close IN;




open IN, $gfa_file || die "fail open $gfa_file";
while (<IN>) {	
	my $line = $_;
	
	if (! exists $DeletedLines{$line}) {
	
		print $line;
	}
}