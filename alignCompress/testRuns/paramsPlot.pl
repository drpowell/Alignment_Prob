#!/usr/bin/perl -w

use strict;

open(F, "> $$.plotdata");
my @keys;

while (<>) {
  next if (/^\#/);

  if (/^PRSS: s1=(\d+) s2=(\d+) parent=(\d+) mutates=(\d+) p=(\S+) s=(\d+)/) {
    next;
  }

  if (/^AlignCompress \(sum=(true|false)\): s1=(\d+) s2=(\d+) parent=(\d+) mutates=(\d+).*?params:(.*)$/) {
    my $mutates = $5;
    my $sum = ($1 eq 'true');
    my $p = $6;
    my $r = {};
    for my $a (split /\s+/,$p) {
      my($k,$v) = split /=/, $a, 2;
      $r->{$k} = $v;
    }

    next if (!$sum);

    if (!@keys) { @keys = sort keys %$r; };

    print F "$mutates ";
    for my $k (sort keys %$r) {
      print F $r->{$k}," ";
    }
    print F "\n";

    next;
  }


  print STDERR "Unknown line: $_";
}
close(F);

my $cmd;

for my $i (2) {
  if ($cmd) {
    $cmd .= ", ";
  } else {
    $cmd = "plot ";
  }
#  $cmd .= sprintf "'$$.plotdata' using 1:%d title '%s'", $i+2, $keys[$i];
  $cmd .= sprintf "'$$.plotdata' thru \$5/\$4 title '2/3'";
}
$cmd .= "\n";

$cmd .= "pause -1\n";
open(F, "| gnuplot -") or (die "Can't run gnuplot");
{ my $oldfh = select F; $|=1; select $oldfh; }
print F $cmd;
<STDIN>;
close(F);
#print $cmd;

unlink("$$.plotdata");
