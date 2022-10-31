#!/usr/bin/perl

=head1 Name

turningpoint_by_lineartransform.pl  -- get the turning point by adjusting  and linear transforming

=head1 Description

detect the turning point automatically for a given set of contact values among contigs

=head1 Version

  Author: Fan Wei, fanwei@caas.cn
  Version: 1.0,  Date: 2021/7/23
  Note:

=head1 Usage
  
  turningpoint_by_lineartransform.pl [options] <CtgContact_file>
  --verbose   output running progress information to screen  
  --help            output help information to screen  

=head1 Example

  perl turningpoint_by_lineartransform.pl formal_100000_iced.matrix.revised.100000.30.CtgContact > formal_100000_iced.matrix.revised.100000.30.CtgContact.adjustTransform 2> formal_100000_iced.matrix.revised.100000.30.CtgContact.turningPoint

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


use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;

my $contact_value_file = shift; 

##matrix for anticlockwise rotate 45 degree
my $A = [
	[sqrt(2)/2, -sqrt(2)/2],
	[sqrt(2)/2, sqrt(2)/2]
];


my @Data;   ## raw data, can be sorted or not sorted
my @AdjustData;  ## Ajusted data divided by a times
my @TransData;  ##linear transformed data by anticlockwise rotate 45 degree
my $Count;
my $Max;
my $Times;

##Load the raw data from input file
my $rank = 1;
open IN, $contact_value_file || die "fail open $contact_value_file\n";
while (<IN>) {
	chomp;
	next if(/^\#/);
	my @t = split /\s+/;
	my $value = $t[2];
	
	my $X = [[$rank],[$value]];
	$Max = $value if($value > $Max);

	push @Data, $X;
	
	$rank ++;
}
close IN;
#print Dumper \@Data;

## Ajusted data divided by a times
$Count = @Data;
$Times = $Max / $Count;
foreach my $p (@Data) {
	#print $p->[1][0],"\n";
	my $order = $p->[0][0];
	my $value = $p->[1][0] / $Times;
	push @AdjustData, [[$order],[$value]];
	#print $p->[0][0],"\t",$p->[1][0],"\n";
}

#print Dumper \@AdjustData;


##linear transformed data by anticlockwise rotate 45 degree
foreach my $p (@AdjustData) {
	my $adujsted = MxM($A, $p);
	push @TransData, $adujsted;
}
#print Dumper \@TransData;

##get the lowest point in the linear transformed data
my $Lowest = $TransData[0][1][0];  ##assign the inital value
my $Lowest_i = 0;
for (my $i = 0; $i < @TransData; $i ++) {
	my $p = $TransData[$i];
	##print $p->[1][0]."\t".$Lowest."\n";
	if ($p->[1][0] < $Lowest) {
		$Lowest = $p->[1][0];
		$Lowest_i = $i;
	}
}

my $turning_point = $Data[$Lowest_i][1][0];
print STDERR "Turning point: $turning_point\n";

##output all the data
print "#Rank\tContacts\tAdjustRank\tAdjustContacts\tTransRank\tTransContacts\n";
for (my $i = 0; $i < @Data; $i ++) {
	print $Data[$i][0][0]."\t".$Data[$i][1][0]."\t".$AdjustData[$i][0][0]."\t".$AdjustData[$i][1][0]."\t".$TransData[$i][0][0]."\t".$TransData[$i][1][0]."\n";
}


##################################################################################

sub MxM {
	my $M1 = shift;
	my $M2 = shift;
	my $M3;

	##用nr,rq表示矩阵的参数个数
	my $M1_n = @$M1;
	my $M1_r = @{$M1->[0]};
	my $M2_r = @$M2;
	my $M2_q = @{$M2->[0]};

	#print "$M1_n\t$M1_r\n$M2_r\t$M2_q\n";

	for (my $i = 0; $i < $M1_n; $i++) {
		for (my $j = 0; $j < $M2_q; $j++) {
			my $value = 0;
			for (my $k = 0; $k < $M1_r; $k ++) {
				$value += $M1->[$i][$k] * $M2->[$k][$j];
			}
			$M3->[$i][$j] = $value;
		}
	}
	
	return $M3;
}


sub dispM{
	my $C = shift;

	foreach my $p (@$C) {
		foreach my $value (@$p) {
			print "$value\t";
		}
		
	}
	print "\n";
}