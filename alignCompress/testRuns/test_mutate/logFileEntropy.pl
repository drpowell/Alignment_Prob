#!/usr/bin/perl -w

use strict;

my $entProg = "./entropy.pl";

my($archetype, $population);

open(F, shift) or die "Can't open";
my $s1;
while (<F>) { $s1.=$_; last if /;$/;};
eval($s1);

$s1 = "";
while (<F>) { $s1.=$_; last if /;$/;};
eval($s1);


my @sum;
my $num=0;
for my $i (0..$#$archetype) {
  my $str = $archetype->[$i]{SEQ};
  my $r = `echo $str | $entProg 2>/dev/null`;
  printf "STR %d:\n", $i;

  while ( $r =~ m/^(\((\d)-order.*?\((\S+) bits\/char\))$/mg ) {
    print "$1\n";
    $sum[$2]+=$3;
    $num++ if ($2 == 0);
  }
}
print "parent Averages\n";
for my $i (0..$#sum) { print "$i: ",$sum[$i]/$num,"\n"; };

@sum=();
$num=0;
for my $i (0..$#$population) {
  my $str = $population->[$i]{SEQ};
  my $r = `echo $str | $entProg 2>/dev/null`;
  printf "STR %d:\n", $i;

  while ( $r =~ m/^(\((\d)-order.*?\((\S+) bits\/char\))$/mg ) {
    print "$1\n";
    $sum[$2]+=$3;
    $num++ if ($2 == 0);
  }
}
print "child Averages\n";
for my $i (0..$#sum) { print "$i: ",$sum[$i]/$num,"\n"; };
