#!/usr/bin/perl -w

use strict;

my $add=2.10480775184078;
my $mult=0.3;

my($archetype, $population);

open(F, shift) or die "Can't open";
my $s1;
while (<F>) { $s1.=$_; last if /;$/;};
#$s1 =~ s/\$\@archetype/\$parents/;
eval($s1);

$s1 = "";
while (<F>) { $s1.=$_; last if /;$/;};
#$s1 =~ s/\$\@population/\$population/;
eval($s1);

my(@prss, @prss2, @al_all, @al_one, @al_sw);

while (<>) {
  next if (/^\#/);

  if (/^PRSS: s1=(\d+) s2=(\d+) parent=(\d+) mutates=(\d+) p=(\S+) s=(\d+)/) {
    push(@prss,{S1=>$1, S2=>$2, PARENT=>$3, MUTATES=>$4, VAL=>(1-$5)});
    push(@prss2,{S1=>$1, S2=>$2, PARENT=>$3, MUTATES=>$4, VAL=>($6)});

    { my($s1,$s2,$v) = ($1,$2,$6);
      my $i;
      for($i=$#al_one; $i>=0; $i--) {
	if ($al_one[$i]{S1} == $s1 && $al_one[$i]{S2} == $s2) { last; }
      }
      die "No find" if ($i<0);
      my $len = ($archetype->[$s1]{S_LEN} + $archetype->[$s1]{SUB_LEN} + $archetype->[$s1]{E_LEN});
      $len += ($population->[$s2]{S_LEN} + $population->[$s2]{SUB_LEN} + $population->[$s2]{E_LEN});
      my $v2 = -$al_one[$i]{VAL};
      my $vv = -($v2 - $add*$len)/$mult;
      print "len=$len\n";
      printf "s1=$s1 s2=$s2 SW_VAL=%f VAL=%f CVAL=%f\n", $v, $v2, $vv;
    }
    next;
  }

  if (/^AlignCompress \(sum=(true|false)\): s1=(\d+) s2=(\d+) parent=(\d+) mutates=(\d+) r=(\S+) \(([\d.]+)\).* ml=\S+ \(([\d.]+)\) dl=\S+ \(([\d.]+)\)/) {
    my $d = {S1=>$2, S2=>$3, PARENT=>$4, MUTATES=>$5, VAL=>-$7};
    if ($1 eq 'true') {
      push(@al_all,$d);
    } else {
      push(@al_one,$d);
    }
    next;
  }

  print STDERR "Unknown line: $_";
}

