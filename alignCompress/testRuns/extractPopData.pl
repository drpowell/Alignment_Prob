#!/usr/bin/perl -w

use strict;

my $prss;
my $blend=0;

if (defined($ARGV[0]) && $ARGV[0] eq 'blend') {
  shift;
  $blend=1;
}

if (defined($ARGV[0]) && (
    $ARGV[0] eq 'prss' ||          # Standard prss p-val
    $ARGV[0] eq 'prss2' ||         # PRSS raw S-W score
    $ARGV[0] eq 'al_all' ||        # Sum over all alignments
    $ARGV[0] eq 'al_one'	   # Optimal alignment only
   )) {
    $prss = shift;
}

(defined($prss)) or die "Usage: $0 <prss|prss2|al_all|al_one|al_sw>";

my(@prss, @prss2, @al_all, @al_one, @al_sw);
my(@parents, @children);

my $sum_or_blend = ($blend ? "blend" : "sum");

while (<>) {
  next if (/^\#/);

  if (/^PRSS: s1=(\d+) s2=(\d+) parent=(\d+) mutates=(\d+) p=(\S+) s=(\d+)/) {
    push(@prss,{S1=>$1, S2=>$2, PARENT=>$3, MUTATES=>$4, VAL=>$5});
    push(@prss2,{S1=>$1, S2=>$2, PARENT=>$3, MUTATES=>$4, VAL=>$6});
    $parents[$3]++;
    $children[$2]++;
    next;
  }

  if (/^AlignCompress   \s  \($sum_or_blend=(true|false)\): \s
      s1=(\d+) \s s2=(\d+) \s parent=(\d+) \s mutates=(\d+) \s
      r=(\S+)    \s \(  ([-\d.]+)     \)   \s
        al=\S+   \s \(   [-\d.]+      \)   \s
      ( ml=(\S+) \s \(  ([-\d.]+)     \)   \s)?
        dl=\S+   \s \(  ([-\d.]+)     \)           /x) {
#    my $d = {S1=>$2, S2=>$3, PARENT=>$4, MUTATES=>$5, VAL=>-$11};
#    my $d = {S1=>$2, S2=>$3, PARENT=>$4, MUTATES=>$5, VAL=>$7};
#    my $d = {S1=>$2, S2=>$3, PARENT=>$4, MUTATES=>$5, VAL=>$7+$10};
    my $d = {S1=>$2, S2=>$3, PARENT=>$4, MUTATES=>$5, VAL=>$6};
#    my $d = {S1=>$2, S2=>$3, PARENT=>$4, MUTATES=>$5, VAL=>$6+$9};
    if ($1 eq 'true') {
      push(@al_all,$d);
    } else {
      push(@al_one,$d);
    }
    $parents[$4]++;
    next;
  }

  print STDERR "Unknown line: $_";
}

my $numParents  = scalar @parents;
my $numChildren = scalar @children;
print STDERR "numParents=$numParents numChildren=$numChildren\n";

my $arr;

($prss eq 'prss')   && ($arr = \@prss);
($prss eq 'prss2')  && ($arr = \@prss2);
($prss eq 'al_all') && ($arr = \@al_all);
($prss eq 'al_one') && ($arr = \@al_one);
($prss eq 'al_sw')  && ($arr = \@al_sw);

#$arr = [ grep { $_->{MUTATES}<80 } @$arr];

$arr = [sort { my $v = $b->{VAL} <=> $a->{VAL};
	       ($prss eq 'prss' ? -$v : $v);
	     } @$arr];

my $errors = 0;
my $correct = 0;
for my $l (@$arr) {
#  print STDERR join " ", (map { "$_ => $l->{$_}" } keys %$l), ($l->{S1} == $l->{PARENT} ? "   #\n" : "\n");
  if ($l->{PARENT} == $l->{S1}) {
    $correct++;
  } else {
    $errors++;
  }
#  print "$correct $errors\n";
  printf "%f %f\n", $correct/$numChildren, $errors/(($numParents-1)*$numChildren);
}



if (0) {
  # Find VAL for x error rate (errors per query.  No. of queries is the number of parent seq)
  my $x = 1;
  my($c,$e)=(0,0);
  for my $i (0..$#$arr) {
    if ($arr->[$i]{PARENT} == $arr->[$i]{S1}) {
      $c++;
    } else {
      my $eRate = $e/$errors;
      if ( $e/$numParents  <= $x && ($e+1)/$numParents>$x ) {
	#      printf STDERR "%s: Rate of %f%% achieves $c true positives (of $correct possible).\n",$ prss, 100*$e/$errors;
	#      printf STDERR "%s: Corresponding cut-off VAL between %g and %g\n", $prss, $arr->[$i-1]{VAL}, $arr->[$i]{VAL};
	
	my $str;
	($prss eq 'prss')   && ($str = 'prss');
	($prss eq 'prss2')  && ($str = 'S-W score');
	($prss eq 'al_all') && ($str = 'Summed alignments');
	($prss eq 'al_one') && ($str = 'Optimal alignment');
	printf STDERR "%s & %d & %d & %g & %g \\\\\n",
	  $str, $c, $correct, $arr->[$i-1]{VAL}, $arr->[$i]{VAL};
      }
      $e++;
    }
  }
}


{
  # Find false/positive rates for a given cutoff
  my $str;
  ($prss eq 'prss')   && ($str = 'prss');
  ($prss eq 'prss2')  && ($str = 'S-W score');
  ($prss eq 'al_all') && ($str = 'Summed alignments');
  ($prss eq 'al_one') && ($str = 'Optimal alignment');

  my $cutoff;
  ($prss eq 'prss')   && ($cutoff = 0.00001);
  ($prss eq 'prss2')  && (exit);
  ($prss eq 'al_all') && ($cutoff = 3);
  ($prss eq 'al_one') && ($cutoff = 0);
  my($c,$e)=(0,0);
  for my $i (0..$#$arr) {
    if (($prss eq 'prss') && ($arr->[$i-1]{VAL}<=$cutoff && $arr->[$i]{VAL}>$cutoff) ||
	($prss ne 'prss') && ($arr->[$i-1]{VAL}>=$cutoff  && $arr->[$i]{VAL}<$cutoff)) {
      printf STDERR "$str & $c & $e & $correct & $errors \\\\ \n";
    }


    if ($arr->[$i]{PARENT} == $arr->[$i]{S1}) {
      $c++;
    } else {
      $e++;
    }

  }
}
