#!/usr/bin/perl -w

my $p = ($ARGV[0] =~ /^1|0$/ ? shift : 0);

while (<>) {
	if (/^PRSS p=(\S+)/ && $p) {
		print "$1\n";
	}
	
	if (/^AlignCompress (\S+)/ && !$p) {
		print "$1\n";
	}
}
