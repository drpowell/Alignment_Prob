#!/usr/bin/perl -w

use strict;

my $prss;

if (defined($ARGV[0]) && (
    $ARGV[0] eq 'prss' ||
    $ARGV[0] eq 'prss2' ||
    $ARGV[0] eq 'al_all' ||
    $ARGV[0] eq 'al_one' ||
    $ARGV[0] eq 'al_sw')) {
    $prss = shift;
}

(defined($prss)) or die "Usage: $0 <prss|prss2|al_all|al_one|al_sw>";

my(@prss, @prss2, @al_all, @al_one, @al_sw);

while (<>) {
  next if (/^\#/);

  if (/^PRSS: s1=(\d+) s2=(\d+) parent=(\d+) mutates=(\d+) p=(\S+) s=(\d+)/) {
    push(@prss,{S1=>$1, S2=>$2, PARENT=>$3, MUTATES=>$4, VAL=>(1-$5)});
    push(@prss2,{S1=>$1, S2=>$2, PARENT=>$3, MUTATES=>$4, VAL=>($6)});
    next;
  }

  if (/^AlignCompress \(sum=(true|false)\): s1=(\d+) s2=(\d+) parent=(\d+) mutates=(\d+) r=(\S+) \(([\d.]+)\).*dl=\S+ \(([\d.]+)\)/) {
#    my $d = {S1=>$2, S2=>$3, PARENT=>$4, MUTATES=>$5, VAL=>-$8};
#    my $d = {S1=>$2, S2=>$3, PARENT=>$4, MUTATES=>$5, VAL=>$7};
    my $d = {S1=>$2, S2=>$3, PARENT=>$4, MUTATES=>$5, VAL=>$6};
    if ($1 eq 'true') {
      push(@al_all,$d);
    } else {
      push(@al_one,$d);
    }
    next;
  }

  if (/^AlignCompress \(SW\): s1=(\d+) s2=(\d+) parent=(\d+) mutates=(\d+) r=(\S+)/) {
    push(@al_sw,{S1=>$1, S2=>$2, PARENT=>$3, MUTATES=>$4, VAL=>$5});
    next;
  }

  print STDERR "Unknown line: $_";
}

my $arr;

($prss eq 'prss')   && ($arr = \@prss);
($prss eq 'prss2')   && ($arr = \@prss2);
($prss eq 'al_all') && ($arr = \@al_all);
($prss eq 'al_one') && ($arr = \@al_one);
($prss eq 'al_sw') && ($arr = \@al_sw);

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
