#!/usr/bin/perl -w

my $prss;

if (defined($ARGV[0]) && (
    $ARGV[0] eq 'prss' ||
    $ARGV[0] eq 'al_all' ||
    $ARGV[0] eq 'al_one')) {
    $prss = shift;
}

(defined($prss)) or die "Usage: $0 <prss|al_all|al_one>";

while (<>) {
	if ($prss eq 'prss' && /^PRSS:\s+mutates=(\d+) p=(\S+)/) {
		print "$1 $2\n";
	}
	
	if ($prss eq 'al_all' && /^AlignCompress \(sum=true\): mutates=(\d+) r=(\S+)/) {
		print "$1 $2\n";
	}
	
	if ($prss eq 'al_one' && /^AlignCompress \(sum=false\): mutates=(\d+) r=(\S+)/) {
		print "$1 $2\n";
	}
}
