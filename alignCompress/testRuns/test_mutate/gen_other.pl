#!/usr/bin/perl -w
#
# Try to generate data that will produce a false negative from prss.
# Make a biased 0th order model.  Sample an *unusual* subseq
# from this (done by negating the model and sampling).
# Output 2 sequences with the high entropy subseq embeded in unrelated 
# low entropy sequences.

use strict;
use Data::Dumper;

my @alpha = qw(a t g c);

my $m = {};
fill($m);
norm($m);
#model_power($m, 3);

my $m2 = deep_copy($m);
negate_model($m2);
norm($m2);

my $subseq = gen_sequence($m2, 40);
#$subseq = "";

print gen_sequence($m, 100), $subseq, gen_sequence($m, 200),"\n";
print gen_sequence($m, 300), $subseq, gen_sequence($m, 50),"\n";

printf("entropy of subseq: under m=%.3f bits/ch   under m2=%.3f bits/ch\n", 
       entropy($m,$subseq), entropy($m2,$subseq));

printf "Entropy: m=%.3f bits m2=%.3f bits\n", model_entropy($m), model_entropy($m2);
print Dumper($m),"\n";
print Dumper($m2),"\n";

sub gen_sequence {		# Generate a sequence from 1st order markov
  my($probs, $len) = @_;

  my $str;

  for (1..$len) {
    my $last = '';      #substr($str, -1);
    my $p = rand();
    my $c;
    for my $cc (keys %{$probs->{$last}}) {
      $p -= $probs->{$last}{$cc};
      if ($p <= 0) {
	$c = $cc;
	last;
      }
    }
    die "No char decided upon" if !defined($c);
    $str .= $c;
  }
  return $str;
}


sub fill {
  my($m) = @_;
  for my $j (@alpha) {
    $m->{''}{$j} = rand();
  }
}

sub norm {
  my($m) = @_;
  my $s = 0;
  for my $k (@alpha) {	
    $s += $m->{''}{$k};
  }
  for my $k (@alpha) {	
    $m->{''}{$k} /= $s;
  }
}

# Calculate the entropy of $str under model $probs (assumes a 0th order model!)
sub entropy {
  my($probs, $str) = @_;

  return undef if (length($str)==0);

  my $s = 0;
  for my $i (0..length($str)-1) {
    my $p += $probs->{''}{substr($str,$i,1)};
    $s += -log2($p);
  }
  return $s/length($str);
}

# Calculate the entropy of a 0th order Markov model.
sub model_entropy {
  my($probs) = @_;
  my $e = 0;
  for my $k1 (keys %$probs) {
    my $e2 = 0;
    for my $k2 (keys %{$probs->{$k1}}) {
      $e2 += -$probs->{$k1}{$k2}*log2($probs->{$k1}{$k2});
    }
#    $e += $e2/4;
    $e += $e2;
  }
  return $e;
}

sub deep_copy {
  my($m) = @_;
  my $r;
  if (ref($m) eq 'HASH') {
    for my $k (keys %$m) { $r->{deep_copy($k)} = deep_copy($m->{$k}) };
  } elsif (ref($m) eq 'ARRAY') {
    for my $i (0..$#$m) { $r->[$i] = deep_copy($m->[$i]); };
  } elsif (ref($m) eq 'SCALAR') {
    $r = $$m;
    $r = \$r;
  } elsif (ref($m) eq '') {
    $r = $m;
  }
  return $r;
}

sub negate_model {
  my($probs) = @_;
  for my $k1 (keys %$probs) {
    for my $k2 (keys %{$probs->{$k1}}) {
      $probs->{$k1}{$k2} = 1-$probs->{$k1}{$k2};
    }
  }
}


# Raise all probabilities in the model to a power, and normalise.
sub model_power {
  my($probs, $pwr) = @_;
  for my $k1 (keys %$probs) {
    my $s;
    for my $k2 (keys %{$probs->{$k1}}) {
      $probs->{$k1}{$k2} **= $pwr;
      $s += $probs->{$k1}{$k2};
    }
    for my $k2 (keys %{$probs->{$k1}}) {
      $probs->{$k1}{$k2} /= $s;
    }
  }
}


sub log2 {
  warn "Can't take log of 0" if ($_[0]==0);
  return log($_[0])/log(2); 
}

sub min {
  my($i) = pop;
  for (@_) { ($_<$i) && ($i=$_); }
  $i;
}

sub max {
  my($i) = pop;
  for (@_) { ($_>$i) && ($i=$_); }
  $i;
}
