#!/usr/bin/perl -w

my $prss;

if (defined($ARGV[0]) && (
    $ARGV[0] eq 'prss' ||
    $ARGV[0] eq 'al_all' ||
    $ARGV[0] eq 'al_one')) {
    $prss = shift;
}

(defined($prss)) or die "Usage: $0 <prss|al_all|al_one>";

my(@prss, @al_all, @al_one);

while (<>) {
	if (/^PRSS:\s+mutates=(-?\d+) p=(\S+)/) {
		push(@prss,[$1,$2]);
	}
	
	if (/^AlignCompress \(sum=true\): mutates=(-?\d+) r=(\S+)/) {
		push(@al_all,[$1,$2]);
	}
	
	if (/^AlignCompress \(sum=false\): mutates=(-?\d+) r=(\S+)/) {
		push(@al_one,[$1,$2]);
	}
}

my $arr;

($prss eq 'prss')   && ($arr = \@prss);
($prss eq 'al_all') && ($arr = \@al_all);
($prss eq 'al_one') && ($arr = \@al_one);

my $last;
my $sum = 0;
my $num = 0;
for my $i (0..$#$arr) {
#	print $arr->[$i][0], " ", $arr->[$i][1], "\n";
#	next;

	if (defined($last) && $last != $arr->[$i][0]) {
		print $last, " ", $sum/$num, "\n";
		$num=0;
		$sum=0;
	}

        $num++;
	$sum += $arr->[$i][1];
	$last = $arr->[$i][0];
}
print $last, " ", $sum/$num, "\n" if ($num);
