## convert GFA < 1.2 format to GFA > 1.2 format
## use Jump to replace Link
use strict;

while (<>) {
	if (/^L\s+/) {
		chomp;
		my @t = split /\s+/;
		$t[0] = "J";
		$t[5] = "*";
		my $line = join("\t",@t);
		print $line."\n";
	}else{
		print;
	}

}