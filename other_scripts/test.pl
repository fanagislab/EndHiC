use strict;


my $str = "ptg000010l-;ptg000025l+";
my $new = rc_ctgs_str($str);

print $str."\n";
print $new."\n";

my $pat = "ptg000010-";

if ($str =~ /$pat/) {
	print "coo\n";
}

sub rc_ctgs_str{
	my $ctgs_str = shift;
	$ctgs_str =~ tr/+-/-+/;
	my @t = split /;/, $ctgs_str;
	@t = reverse @t;
	$ctgs_str = join(";", @t);
	return $ctgs_str;
}

