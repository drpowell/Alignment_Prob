#!/usr/bin/perl -w
#
# Count the occurances of numbers is $col.  If floats then bin them to $acc

my $col = shift;

my $acc = 1;
if (@ARGV>1) {
  $acc = shift;
  print STDERR "acc=$acc\n";
}

my %c;

while (<>) {
  my $d = (split /\s+/)[$col];
  $d = $d*$acc;
  $d += ($d>0 ? 0.5 : -0.5);
  $d = int($d)/$acc;
  $c{$d}++;
}

for my $k (sort { $a<=>$b } keys %c) {
  print $k," ",$c{$k},"\n";
}
