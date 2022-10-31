#!/usr/bin/perl

=head1 Name

shortCtgs_to_cluster.pl  -- assign the short contigs into cluster 

=head1 Description

This program consider the max contact as well as the max / second times to cluster, in order to assign the short contigs to clusters

=head1 Version

  Author: Fan Wei, fanwei@caas.cn
  Version: 1.0,  Date: 2021/7/6
  Note:

=head1 Usage
  
  shortCtgs_to_cluster.pl [options] <cluster_file> <contig_length_file> <ctgContact_file>
  --contact <float>  cutoff for the contact value of the max cluster, default=0
  --times <float>  cutoff for the max contact / second contact, default=1
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

  perl ../../shortCtgs_to_cluster.pl --contact 500 --times 2  formal_1000000_iced.matrix.ctgcontact.gfa.cluster  contigs_200kb.len  formal_1000000_iced.matrix.ctgcontact > contigs_200kb.len.cluster
  

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;


##get options from command line into variables and set default values
my $max_contact_cutoff;
my $times_max_second;
my ($Verbose,$Help);
GetOptions(
	"contact:f"=>\$max_contact_cutoff,
	"times:f"=>\$times_max_second,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$max_contact_cutoff ||= 0;
$times_max_second ||= 1;
die `pod2text $0` if (@ARGV == 0 || $Help);

#print STDERR "max_contact_cutoff:  $max_contact_cutoff\n";
#print STDERR "max / second cutoff: $times_max_second\n";

my $cluster_result_file = shift;
my $contig_len_file = shift;
my $ctgContact_file = shift;  ##from ctgEnd or maxBin, either is OK, only keep the largest contact for each contig pair

my @Cluster;
my @ClusterId;
my %ClusteredCtgs;
my @Contig;
my %Contact;
my %CtgLen;


open IN, $cluster_result_file || die "fail open $cluster_result_file";
while (<IN>) {
	next if(/^\#/);
	chomp;
	my @t = split /\s+/;
	my $cluster_id = $t[0];
	my $ctg_count = $t[1];
	my $cluster_len = $t[2];
	my $robustness = $t[3];
	my $ctg_str = $t[4];
	
	my @ary;
	while ($ctg_str =~ /(\w+)[+-]/g) {
		push @ary, $1;
		$ClusteredCtgs{$1} = 1;
	}
	
	push @Cluster, \@ary;
	push @ClusterId, $cluster_id;

}
close IN;

#print Dumper \@ClusterId;
#print Dumper \%ClusteredCtgs;


open IN, $contig_len_file || die "fail open $contig_len_file\n";
while (<IN>) {
	chomp;
	my @t = split /\s+/;
	my $ctg_id = $t[0];
	my $ctg_len = $t[1];
	$CtgLen{$ctg_id} = $ctg_len;
	push @Contig, $ctg_id if(! exists $ClusteredCtgs{$ctg_id});
}
close IN;

#print Dumper \@Contig;



open IN, $ctgContact_file || die "fail open $ctgContact_file\n";
while (<IN>) {
	next if(/^\#/);
	chomp;
	my @t = split /\s+/;
	my $ctg1 = $t[0];
	my $ctg2 = $t[1];
	my $contact = $t[2];  #normalized contact
	
	##only keep the max value from head vs head, head vs tail, tail vs head, tail vs tail
	if (!exists $Contact{$ctg1}{$ctg2} || $Contact{$ctg1}{$ctg2} < $contact) {
		$Contact{$ctg1}{$ctg2} = $contact;
	}

	if (!exists $Contact{$ctg2}{$ctg1} || $Contact{$ctg2}{$ctg1} < $contact) {
		$Contact{$ctg2}{$ctg1} = $contact;
	}
	

}
close IN;

#print Dumper \%Contact;

print "#ContigId\tContigLen\tClusterId\tContact\tTimes\tmax_cluster_id:contig_id:contact_value\tsecond_cluster_id:contig_id:contact_value\n";
my @Output;
foreach my $this_ctg (@Contig) {
	
	my @stored;
	
	for (my $i = 0; $i < @Cluster; $i ++) {
		my $p = $Cluster[$i];
		my $max_contig = "";
		my $max_contact = 0;
		foreach my $cluter_ctg (@$p) {
			#print $this_ctg."\t".$cluter_ctg."\n";
			if(exists $Contact{$this_ctg}{$cluter_ctg}){
				my $contact_value = $Contact{$this_ctg}{$cluter_ctg};
				if ($contact_value > $max_contact) {
					$max_contig = $cluter_ctg;
					$max_contact = $contact_value;
				}
			}
		}

		push @stored, [$max_contact, $max_contig, $i];
	}

	@stored = sort {$b->[0] <=> $a->[0]} @stored;
	
	my ($max_contact , $max_contig , $max_cluster ) =  ($stored[0][0],   $stored[0][1],  $stored[0][2]);
	my ($second_contact , $second_contig , $second_cluster ) =  ($stored[1][0],   $stored[1][1],  $stored[1][2]);
	my $max_cluster_id = $ClusterId[$max_cluster];
	my $second_cluster_id = $ClusterId[$second_cluster];
		
	my $assigned_cluster = "Unclustered";
	$max_contact = 1 if($max_contact == 0);
	$second_contact = $max_contact if($second_contact == 0);
	my $times = sprintf("%.2f", $max_contact / $second_contact);
	if ($max_contact >= $max_contact_cutoff && $times >= $times_max_second) {
		$assigned_cluster = $max_cluster_id;
	}

	
	push @Output, [$this_ctg, $CtgLen{$this_ctg}, $assigned_cluster, $max_contact, $times, "$max_cluster_id:$max_contig:$max_contact", "$second_cluster_id:$second_contig:$second_contact"];
}

##output the sorted results
@Output = sort {$a->[2] cmp $b->[2]} @Output;
foreach my $p (@Output) {
	my $line = join("\t", @$p);
	print $line."\n";
}



