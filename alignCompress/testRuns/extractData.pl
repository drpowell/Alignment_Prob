#!/usr/bin/perl -w

my $prss;

if (defined($ARGV[0]) && (
    $ARGV[0] eq 'prss' ||
    $ARGV[0] eq 'al_all' ||
    $ARGV[0] eq 'al_one' ||
    $ARGV[0] eq 'al_all_blend' ||
    $ARGV[0] eq 'al_one_blend')) {
    $prss = shift;
}

(defined($prss)) or die "Usage: $0 <prss|al_all|al_one|al_all_blend|al_one_blend>";

my(@prss, @al_all, @al_one, @al_all_blend, @al_one_blend);

while (<>) {
	next if (/^#/);
	if (/^PRSS:\s+mutates=(-?\d+) p=(\S+)/) {
		push(@prss,[$1,$2]);
		next;
	}
	
	if (/^AlignCompress \((?:blend=false )?sum=true\): mutates=(-?\d+) r=(\S+)/) {
		push(@al_all,[$1,$2]);
		next;
	}
	
	if (/^AlignCompress \((?:blend=false )?sum=false\): mutates=(-?\d+) r=(\S+)/) {
		push(@al_one,[$1,$2]);
		next;
	}

	if (/^AlignCompress \(blend=true sum=true\): mutates=(-?\d+) r=(\S+)/) {
		push(@al_all_blend,[$1,$2]);
		next;
	}
	
	if (/^AlignCompress \(blend=true sum=false\): mutates=(-?\d+) r=(\S+)/) {
		push(@al_one_blend,[$1,$2]);
		next;
	}
	print STDERR "Unknown line: $_";
}

my $arr;

($prss eq 'prss')   && ($arr = \@prss);
($prss eq 'al_all') && ($arr = \@al_all);
($prss eq 'al_one') && ($arr = \@al_one);
($prss eq 'al_all_blend') && ($arr = \@al_all_blend);
($prss eq 'al_one_blend') && ($arr = \@al_one_blend);

my $last;
my $sum = 0;
my $num = 0;
my($low,$high);
for my $i (0..$#$arr) {
#	print $arr->[$i][0], " ", $arr->[$i][1], "\n";
#	next;

	if (defined($last) && $last != $arr->[$i][0]) {
		print $last, " ", $sum/$num, " ", $low, " ", $high, "\n";
		$num=0;
		$sum=0;
	}

	$low  = $arr->[$i][1] if ($num==0 || $arr->[$i][1]<$low);
	$high = $arr->[$i][1] if ($num==0 || $arr->[$i][1]>$high);

	$sum += $arr->[$i][1];
	$last = $arr->[$i][0];

        $num++;
}
print $last, " ", $sum/$num, " ", $low, " ", $high, "\n" if ($num);
