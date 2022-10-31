#!/usr/bin/perl

=head1 Name

summarize_endhic_results.pl  -- summarize and analyze a set of EndHiC results 

=head1 Description

Perform robustness analysis for all cluster types and find out stable cluster types

=head1 Version

  Author: Fan Wei, fanwei@caas.cn
  Version: 1.0,  Date: 2021/7/23
  Note:

=head1 Usage
  
  summarize_endhic_results.pl [options] <*.reciprocalMax.gfa.cluster.order.orient>
  --ctglenfile <str>  contig id and length file
  --ctglencut <int>   length cutoff for used large contigs, default=0
  --binsize <int>   bin size used in the hicpro result files, default=100000
  --clustermark <str>    mark character for cluster ID, default=A
  --help            output help information to screen  

=head1 Example

  perl summerize_endhic_results.pl  --binsize 100000  *.reciprocalMax.gfa.cluster.order.orient  > z.EndHiC.results.summary.and.analysist

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;


##get options from command line into variables and set default values
my $Contig_used_file;
my $Contig_len_cutoff;
my $Times;
my $BinSize;
my $ClusterMark;
my ($Verbose,$Help);
GetOptions(
	"ctglenfile:s"=>\$Contig_used_file,
	"ctglencut:i"=>\$Contig_len_cutoff,
	"binsize:i"=>\$BinSize,
	"clustermark:s"=>\$ClusterMark,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$BinSize ||= 100000;
$Contig_len_cutoff ||= 0;
$ClusterMark ||= "";
die `pod2text $0` if (@ARGV == 0 || $Help);

my @InputFiles = @ARGV;

my %ClusterData;

my %Summary;

my %UsedContigs;

##get contig id and length data
open IN, $Contig_used_file || die "fail open $Contig_used_file\n";
while (<IN>) {
	if(/(\S+)\s+(\d+)/){
		$UsedContigs{$1} = $2;
	}
}
close IN;


my $FileNum = @InputFiles;

foreach my $file (@InputFiles) {
	my $cluster_num = 0;
	open IN, $file || die "fail open $file";
	while (<IN>) {
		next if(/^\#/ || /Unclustered/);
		$cluster_num ++;
		chomp;
		my ($cluster_id, $ctg_num, $cluster_len, $robustness, $ctg_str) = split /\s+/;
	    
		##the $ctg_str already has order, so no need to consider its reverse complmentary strand form
		$ClusterData{$ctg_str} ++;   
	}
	close IN;
	
	#print "In $file\nCluster number: $cluster_num\n\n";
	##formal_100000.matrix.100000.5.CtgContact.overCutoff.2.5.reciprocalMax.gfa.cluster.order.orient
	my $binNum = $1 if($file =~ /\.(\d+)\.CtgContact.overCutoff/ || $file =~ /\.(\w+)\.overCutoff/);
	my $times = $1 if($file =~ /overCutoff\.(\S+)\.reciprocalMax.gfa/);
	my $rawIced = ($file =~ /_iced\.matrix/) ? "Iced" : "Raw";
	$Summary{$binNum}{$times}{$rawIced} = $cluster_num;
}

print "Number of clusters under each condition:\n";

foreach my $binNum (sort {$a <=> $b} keys %Summary) {
	my $binNum_p = $Summary{$binNum};
	print "\nContig end size: $BinSize x $binNum\n";
	
	my $times_str = "Times:";
	my $raw_str =   "Raw:  ";
	my $iced_str =  "Iced: ";
	foreach my $times (sort keys %$binNum_p) {
		my $times_p = $binNum_p->{$times};
		my $raw = $times_p->{"Raw"};
		my $iced = $times_p->{"Iced"};
		$times_str .= "\t".$times;
		$raw_str .= "\t".$raw;
		$iced_str .= "\t".$iced;
	}
	print $times_str."\n";
	print $raw_str."\n";
	print $iced_str."\n";
}

print "\n------------------------------------------------------------------------\n";


my $allcluster_number = keys %ClusterData;
print "\nNumber of all Cluster units: $allcluster_number\n\n";

print "#ctgs_order_orientation\trobustness[max:$FileNum]\tcluster_length\n";


foreach my $cluster  (sort {$ClusterData{$b} <=> $ClusterData{$a}} keys %ClusterData) {
	my $frequency = $ClusterData{$cluster};
	my $cluster_len = cal_cluster_len($cluster);
	print "$cluster\t$frequency\t$cluster_len\n";
}

print "\n------------------------------------------------------------------------\n";


my %CombineCluster;
my %CombineClusterTypes;

foreach my $cluster  (sort {$ClusterData{$b} <=> $ClusterData{$a}} keys %ClusterData) {
	
	##ptg000055l-;ptg000033l-;ptg000058l-
	my $frequency = $ClusterData{$cluster};
	my @ctgs = split /;/, $cluster;
	my @filtered_ctgs;
	foreach my $ctg_str (@ctgs) {
		my $ctg_id = $ctg_str;
		$ctg_id =~ s/[+-]$//;
		my $ctg_len = $UsedContigs{$ctg_id};
		push @filtered_ctgs, $ctg_str if($ctg_len >= $Contig_len_cutoff);
	}
	my $filtered_ctgs_cluster = join(";",@filtered_ctgs);
	
	if (@filtered_ctgs == 1) {
		if ($filtered_ctgs_cluster =~ /-/) {
			$filtered_ctgs_cluster =~ s/-/+/;
		}
	}
	if (@filtered_ctgs > 1) {
		if ($filtered_ctgs[0] gt $filtered_ctgs[-1]) {
			$filtered_ctgs_cluster = rc_ctgs_str($filtered_ctgs_cluster);
		}
	}

	if (@filtered_ctgs >= 1) {
		$CombineCluster{$filtered_ctgs_cluster} += $frequency;
		push @{$CombineClusterTypes{$filtered_ctgs_cluster}}, $cluster, $frequency;
	}
}

my $filtered_cluster_number = keys %CombineCluster;
print "\nNumber of combined Cluster units: $filtered_cluster_number\n\n";

print "#ctgs_order_orientation\trobustness[max:$FileNum]\tcluster_length\tcombined_cluster_number\tcombined_cluster_details\n";

my %AppearedCtgs;
my @Confident_cluster; ##for high robustness to low robustness
foreach my $cluster  (sort {$CombineCluster{$b} <=> $CombineCluster{$a}} keys %CombineCluster) {
	my $frequency = $CombineCluster{$cluster};
	my $cluster_types_p = $CombineClusterTypes{$cluster};
	my $cluster_types = @$cluster_types_p / 2;
	my $cluster_types_str = join(" ",@$cluster_types_p); 
	$cluster_types_str = "[$cluster_types_str]";
	

	my $cluster_len = cal_cluster_len($cluster);
	
	my @ctgs = split /;/, $cluster;
	
	my $is_existed = 0;
	foreach my $ctg (@ctgs) {
		$ctg =~ s/[+-]//;
		if (exists $AppearedCtgs{$ctg}) {
			$is_existed = 1;
		}
	}
	if ($is_existed == 0) {  ##non-redundant contigs set
		push @Confident_cluster, [$cluster_len, $frequency,$cluster, $cluster_types_str];
		foreach my $ctg (@ctgs) {
			$ctg =~ s/[+-]//;
			$AppearedCtgs{$ctg} = 1;
		}
	}		
	
	print "$cluster\t$frequency\t$cluster_len\t$cluster_types\t$cluster_types_str\n" ;

}

print "\n------------------------------------------------------------------------\n";


my $confidentcluster_number = @Confident_cluster;
print "\nNumber of stable cluster units: $confidentcluster_number\n\n";

print "#Cluster_id\tcontig_count\tcluster_length\trobustness[max:$FileNum]\tcontigs_order_orientation\trobustness_for_this_cluster_type[max:$FileNum]\n";
#print STDERR "#Cluster_id\tcontig_count\tcluster_length\trobustness[max:$FileNum]\tcontigs_order_orientation\n";


my %Confident_contigs;
my %Confident_contigs_appearance;
my @Confident_cluster_new;
foreach my $p (sort {$b->[0] <=> $a->[0]} @Confident_cluster) { #
	my ($cluster_len, $frequency,$cluster, $cluster_types_str) = @$p;
	$cluster_types_str =~ s/\[//;
	$cluster_types_str =~ s/\]//;
	my @t = split /\s+/, $cluster_types_str;
	
	##add back the small contigs, and get the longest cluster including as max number of contigs as possible
	my $max_cluster = $t[0];
	my $max_cluster_freq = $t[1];
	for (my $i=2; $i<@t; $i+=2) {
		my $cluster = $t[$i];
		my $cluster_freq = $t[$i+1];
		if (length($cluster) > length($max_cluster)) {
			$max_cluster = $cluster;
			$max_cluster_freq = $cluster_freq;
		}
	}
	

	$cluster_len = cal_cluster_len($max_cluster);
	my @ctgs = split /;/, $max_cluster;
	foreach my $ctgId (@ctgs) {
		$ctgId =~ s/[+-]//;
		$Confident_contigs{$ctgId} = $UsedContigs{$ctgId};
		$Confident_contigs_appearance{$ctgId} ++;
	}
	my $ctgs_num = @ctgs;
	push @Confident_cluster_new, [$ctgs_num,$cluster_len,$frequency,$max_cluster, $max_cluster_freq];
	
	
}

my @Confident_cluster_final; ##without short redundant contigs
my $loop_id = "01";
foreach my $p (sort {$b->[1] <=> $a->[1]} @Confident_cluster_new) {
	
	push @Confident_cluster_final, ["Cluster_$ClusterMark$loop_id", @$p];

	my $line = join("\t", @$p);
	print "Cluster_$ClusterMark$loop_id\t$line\n";
	#print STDERR "Cluster_$ClusterMark$loop_id\t$line\n";
	$loop_id ++;
}


print "\n------------------------------------------------------------------------\n";


print "\nNumber of stable cluster units (no redundant contigs): $confidentcluster_number\n\n";

print "#Cluster_id\tcontig_count\tcluster_length\trobustness[max:$FileNum]\tcontigs_order_orientation\n";
print STDERR "#Cluster_id\tcontig_count\tcluster_length\trobustness[max:$FileNum]\tcontigs_order_orientation\n";

my %RedundantList;
foreach my $ctgId (sort keys %Confident_contigs_appearance) {
	my $appearance = $Confident_contigs_appearance{$ctgId};
	$RedundantList{$ctgId} = 1 if($appearance > 1);
}

#print STDERR Dumper \%RedundantList;


##remove the redundant short contigs from each stable clusters
foreach my $ctgId (sort keys %RedundantList) {
	
	my $max_freq = 0;
	my $max_freq_id = 0;
	for (my $i=0; $i<@Confident_cluster_final; $i++) {
		my $ctgs_str = $Confident_cluster_final[$i][4];
		next if($ctgs_str !~ /$ctgId[+-]/);
		my $cluster_type_freq = $Confident_cluster_final[$i][5];
		if ($max_freq < $cluster_type_freq) {
			$max_freq = $cluster_type_freq;
			$max_freq_id = $i;
		}
	}
	#print STDERR "$max_freq\t$max_freq_id\n";
	for (my $i=0; $i<@Confident_cluster_final; $i++) {
		my $ctgs_str = $Confident_cluster_final[$i][4];
		next if($ctgs_str !~ /$ctgId[+-]/);
		if ($i != $max_freq_id) {
			$Confident_cluster_final[$i][2] -= $UsedContigs{$ctgId};
			$Confident_cluster_final[$i][4] =~ s/$ctgId[+-];?//;
		}
	}
}

my $loop_id = "01";
foreach my $p (sort {$b->[2] <=> $a->[2]} @Confident_cluster_final) {
	$p->[0] = "Cluster_$ClusterMark$loop_id";
	pop @$p;
	my $line = join("\t",@$p);
	$line =~ s/;$//;
	print $line."\n";
	print STDERR $line."\n";
	$loop_id ++;
}

#print STDERR Dumper \@Confident_cluster_final;



print "\n------------------------------------------------------------------------\n";


print "\nIncluding contigs:\n\n";

print "#ctgId\tctg_len\tctg_appearance\tnote\n";
my $total_contig_num = 0;
my $total_contig_len = 0;
foreach my $ctgId (sort {$Confident_contigs{$b} <=> $Confident_contigs{$a}} keys %Confident_contigs) {
	my $ctgId_len = $Confident_contigs{$ctgId};
	my $ctgId_appearance = $Confident_contigs_appearance{$ctgId};
	$ctgId_appearance .= "\tNote: this contig appears more than one times in clusters, redundance have been removed" if($ctgId_appearance > 1);
	print "$ctgId\t$ctgId_len\t$ctgId_appearance\n";
	$total_contig_num ++;
	$total_contig_len += $ctgId_len;
}

print "\nTotal Contig number: $total_contig_num\n";
print "Total contig length: $total_contig_len";

###########################################################################


sub rc_ctgs_str{
	my $ctgs_str = shift;
	$ctgs_str =~ tr/+-/-+/;
	my @t = split /;/, $ctgs_str;
	@t = reverse @t;
	$ctgs_str = join(";", @t);
	return $ctgs_str;
}


sub cal_cluster_len{
	my $cluster = shift;
	my $cluster_len = 0;
	my @ctgs = split /;/, $cluster;
	foreach my $ctgId (@ctgs) {
		$ctgId =~ s/[+-]//;
		$cluster_len += $UsedContigs{$ctgId};
	}
	return $cluster_len;
}

