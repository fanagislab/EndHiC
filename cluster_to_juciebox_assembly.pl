#!/usr/bin/perl -w
use strict;
use Data::Dumper;

if(!@ARGV) {
	print "\nThis program can produce .assmelby file for Juicebox visualization, based on a file of contig len and EndHiC cluster.\n";
	print "Usage:\n\tperl cluster_to_juicebox_assembly.pl contigs.fa.len z.EndHiC.A.results.summary.cluster > draft.assembly\n\n";
	exit(1);
}

my $len_f = shift;
my $cls_f = shift;

my %Len;
my $cnt = 0;
my %Id;
my %Id2;
open IN, $len_f || die "fail open $len_f\n";
while(<IN>) {
	chomp;

	my ($id, $len) = (split /\t/,$_)[0, 1];
	$cnt++;
	$Len{$cnt} = $len;
	$Id{$cnt} = $id;
	$Id2{$id} = $cnt;
}
close IN;
#print Dumper \%Len;
#print Dumper \%Id;

foreach my $id (sort {$a<=>$b}keys %Len) {
	my $len = $Len{$id};
	my $ctg = $Id{$id};
	print ">$ctg $id $len\n";
}

++$cnt;
print ">hic_gap_$cnt $cnt 1000\n";

my %Anchored;
open IN1, $cls_f || die "fail open $cls_f\n";
while(<IN1>) {
	chomp;

	next if($_ =~ /^#/);
	my $ctg = (split /\t/,$_)[4];
	my @contigs = split /;/,$ctg;
	my $out = "";
	foreach my $c (@contigs) {
		my $id = $c;
		$id =~ s/[+-]//;
		$Anchored{$id} = 1;

		my $dire = $1 if($c =~ /(\S)$/);
		my $num = $Id2{$id};		# the number cooresponsible to the contig
		$out .= " $dire$num $cnt";
	}
	$out =~ s/^\s//;
	$out =~ s/\s$cnt$//;
	$out =~ s/\+//g;
	print "$out\n";
}
close IN1;
#print Dumper \%Anchored;

foreach my $ctg (sort keys %Id2) {
	if(!exists $Anchored{$ctg}) {
		print "$Id2{$ctg}\n";
	}
}
