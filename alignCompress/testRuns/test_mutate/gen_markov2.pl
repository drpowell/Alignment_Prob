#!/usr/bin/perl -w

# Generate a 2nd order Markov Model such that the 1st order statistics are uniform.

use strict;
use Data::Dumper;

my @alpha = qw(a t g c);

my $m = {};
fill($m);
while (norm($m) > 1E-10) { };
model_power($m, 10);
while (norm($m) > 1E-10) { };


print gen_sequence($m,1000),"\n";
print gen_sequence($m,1000),"\n";
print Dumper($m),"\n";

sub gen_sequence {		# Generate a sequence from 2st order markov
  my($probs, $len) = @_;

  my $str; # Pick any 2 starting chars
  $str .= $alpha[rand(@alpha)] . $alpha[rand(@alpha)];

  for (3..$len) {
    my $last = substr($str, -2);
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
  for my $i (@alpha) {
    for my $j (@alpha) {
      for my $k (@alpha) {
	$m->{$i.$j}{$k} = rand();
      }
    }
  }
}

sub norm {
  my($m) = @_;
  my($max) = 0;
  for my $i (@alpha) {
    for my $j (@alpha) {
      my $s = 0;
      for my $k (@alpha) {	
	$s += $m->{$i.$j}{$k};
      }
      for my $k (@alpha) {	
	$m->{$i.$j}{$k} /= $s;
      }
      $max = max($max, abs($s-1));
    }
  }

  for my $j (@alpha) {
    for my $k (@alpha) {
      my $s = 0;
      for my $i (@alpha) {	
	$s += $m->{$i.$j}{$k};
      }
      for my $i (@alpha) {	
	$m->{$i.$j}{$k} /= $s;
      }
      $max = max($max, abs($s-1));
    }
  }

  return $max;
}

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


sub max {
  my($i);
  for my $v (@_) {
    (!defined($i) || $v>$i) && ($i=$v);
  }
  return $i;
}
