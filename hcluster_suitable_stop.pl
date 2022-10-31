#!/usr/bin/perl

=head1 Name

hcluster_suitable_stop.pl  -- stop clustering at suitable loop

=head1 Description



=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage
  
  hcluster_suitable_stop.pl <hcluster_result_file>
  --chr_num <int>    expected chromosome number
  --chr_len <int>    minimum chromosome length
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

my $chromosome_number;
my $chromosome_length;
my ($Verbose,$Help);
GetOptions(
	"chr_num:i"=>\$chromosome_number,
	"chr_len:i"=>\$chromosome_length,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$chromosome_number ||= 18;
$chromosome_length ||= 10000000;
die `pod2text $0` if (@ARGV == 0 || $Help);

#print STDERR "Parameters: $chromosome_number\t$chromosome_length\t$anchor_rate\n";

my $hcluster_log_file = shift;
my @DATA;

$/ = "Clustering loop";
open IN, $hcluster_log_file || die "fail $hcluster_log_file";
<IN>;
while (<IN>) {
	chomp;
	s/Clustering finished\s+//;
	#print "\n--------------------------------------------\n";
	my @lines = split /\n+/;
	my $loop_line = shift @lines;
	my $loop_id = $1 if($loop_line =~ /^\s+(\d+)/);
	my $group_num_line = shift @lines;
	#print Dumper \@lines;
	#print $loop_id."\n";
	
	
	my @data;
	foreach  (@lines) {
		my @t = split /\s+/;
		my ($cluster_id, $contig_num,$cluster_len, $robustness, $contigs_str) = @t;
		push @data, [$cluster_id, $contig_num,$cluster_len, $robustness, $contigs_str];
	}
	push @DATA, \@data;

}
close IN;

@DATA = reverse @DATA;


#print Dumper \@DATA;


for (my $i = $chromosome_number - 1; $i < @DATA; $i++) {
	
	my $cluster_number = 0;
	
	my $p = $DATA[$i];
	foreach my $pp (@$p) {
		if($pp->[2] > $chromosome_length){
			$cluster_number ++;
		}
	}
	
	#print STDERR "$i\t$cluster_number\n";
	if ($cluster_number == $chromosome_number) {
		my $loop_id = $i + 1;
		print "#Cluster_id\tcontig_count\tcluster_length\trobustness\tincluded_contigs(no_order_and_orient) [Loop_id: $loop_id]\n";
		foreach my $pp (@$p) {
			my $line = join("\t",@$pp);
			print $line."\n";
		}
		last;
	}

}



####################################################
################### Sub Routines ###################
####################################################
