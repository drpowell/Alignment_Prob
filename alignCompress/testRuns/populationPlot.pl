#!/usr/bin/perl -w

use strict;

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
  next if (/^\#/);

  if (/^PRSS: s1=(\d+) s2=(\d+) parent=(\d+) mutates=(\d+) p=(\S+)/) {
    push(@prss,{S1=>$1, S2=>$2, PARENT=>$3, MUTATES=>$4, VAL=>(1-$5)});
    next;
  }

  if (/^AlignCompress \(sum=true\): s1=(\d+) s2=(\d+) parent=(\d+) mutates=(\d+) r=(\S+)/) {
    push(@al_all,{S1=>$1, S2=>$2, PARENT=>$3, MUTATES=>$4, VAL=>$5});
    next;
  }

  if (/^AlignCompress \(sum=false\): s1=(\d+) s2=(\d+) parent=(\d+) mutates=(\d+) r=(\S+)/) {
    push(@al_one,{S1=>$1, S2=>$2, PARENT=>$3, MUTATES=>$4, VAL=>$5});
    next;
  }
}

my $arr;

($prss eq 'prss')   && ($arr = \@prss);
($prss eq 'al_all') && ($arr = \@al_all);
($prss eq 'al_one') && ($arr = \@al_one);

#$arr = [ grep { $_->{MUTATES}>80 } @$arr];

$arr = [sort { $b->{VAL} <=> $a->{VAL} } @$arr];

my $errors = 0;
my $correct = 0;
for my $l (@$arr) {
#  print join " ", (map { "$_ => $l->{$_}" } keys %$l), "\n";
  if ($l->{PARENT} == $l->{S1}) {
    $correct++;
  } else {
    $errors++;
  }
  print "$correct $errors\n";
}