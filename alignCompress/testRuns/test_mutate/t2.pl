#!/usr/bin/perl -w


# This script attempts to mutate a sequence generated from a 1st order
# Markov model while still keeping the mutated sequence _in_ the model.
# This approach is inspired by discussion with Lloyd and Wallace's proof
# for the Metropolis algorithm (http://www.csse.monash.edu.au/~lloyd/tildeMML/Local/1995-Metropolis/)
#
# Works as follows:
#    Generate sequence s1 using MM.
#
# Similar to t.pl, but has a parameter to control the ratio
# of changes to indels.


use strict;
use Data::Dumper;

use UKK;

my($Pchange) = (0.8);

use vars qw($main_loop);	# Global for main loop counter. So can access in signal handler
$SIG{HUP} = sub {print STDERR "Upto: $main_loop\n";};


my($len_lower, $len_upper) = (100,200);

my $probs = { 'a' => {'a' => 0.11, 't'=>0.09, 'g'=>0.3,  'c'=>0.5},
	      't' => {'a' => 0.09, 't'=>0.11, 'g'=>0.5,  'c'=>0.3},
	      'g' => {'a' => 0.5,  't'=>0.3,  'g'=>0.09, 'c'=>0.11},
	      'c' => {'a' => 0.5,  't'=>0.11, 'g'=>0.09, 'c'=>0.3}
	    };

model_power(0.001);
printf "Model = %s\n", Dumper($probs);
printf "Model entropy = %f\n", model_entropy();

my(%mutate_counts);

#for (1..100000) {
#  print rand_gauss(5,2),"\n";
#}
#exit;

my(@e1,@e2);
my(@l1,@l2);
my(@ed);
for $main_loop (1..1000) {
  my $str = gen_sequence(int(rand($len_upper-$len_lower+1))+$len_lower);
  push(@e1, entropy($str));
  push(@l1, length($str));

#  my $str2 = gen_sequence(int(rand($len_upper-$len_lower+1))+$len_lower);
  my $str2 = mutate($str, 4000);
  push(@e2, entropy($str2));
  push(@l2, length($str2));

  my $ukk = new UKK($str,$str2);
  push(@ed, $ukk->edit_dist);

#  print "$str\n$str2\n";
#  print "@ed\n";
#  exit;
}

printf "Model entropy = %f\n", model_entropy();
printf "sequence1 entropy mean=%f std=%f\n", fit_gauss(\@e1);
printf "sequence1 lengths mean=%f std=%f\n", fit_gauss(\@l1);
printf "sequence2 entropy mean=%f std=%f\n", fit_gauss(\@e2);
printf "sequence2 lengths mean=%f std=%f\n", fit_gauss(\@l2);
printf "edit dist mean=%f std=%f min=%f max=%f\n", fit_gauss(\@ed), min(@ed), max(@ed) if (@ed);
for my $op (sort keys %mutate_counts) { print "$op: $mutate_counts{$op} "; } print "\n";


sub mutate_entropy {
  # This is for 1st order MMs only!
  my($entropy, $str, $len, $op) = @_;

  if ($op->[0] eq 'c') {
    my($p, $newc) = ($op->[1], $op->[2]);
    my($pc, $cc, $nc) = (substr($$str, $p-1, 1),
			 substr($$str, $p,   1),
			 substr($$str, $p+1, 1));
    my($e1,$e2) = (0,0);
    $e1 += -log2($probs->{$pc}{$cc}) if ($p>0);
    $e1 += -log2($probs->{$cc}{$nc}) if ($p<$len-1);
    $e2 += -log2($probs->{$pc}{$newc}) if ($p>0);
    $e2 += -log2($probs->{$newc}{$nc}) if ($p<$len-1);
    return $entropy - $e1 + $e2;
  }

  if ($op->[0] eq 'd') {
    my($p) = ($op->[1]);
    my($pc, $cc, $nc) = (substr($$str, $p-1, 1),
			 substr($$str, $p,   1),
			 substr($$str, $p+1, 1));
    my($e1,$e2) = (0,0);
    $e1 += -log2($probs->{$pc}{$cc}) if ($p>0);
    $e1 += -log2($probs->{$cc}{$nc}) if ($p<$len-1);
    $e2 += -log2($probs->{$pc}{$nc}) if ($p>0 && $p<$len-1);
    return $entropy - $e1 + $e2;
  }

  if ($op->[0] eq 'i') {
    my($p, $newc) = ($op->[1], $op->[2]);
    my($pc, $nc) = (substr($$str, $p-1, 1),
		    substr($$str, $p,   1));
    my($e1,$e2) = (0,0);
    $e1 += -log2($probs->{$pc}{$nc}) if ($p>0 && $p<$len-1);
    $e2 += -log2($probs->{$pc}{$newc}) if ($p>0);
    $e2 += -log2($probs->{$newc}{$nc}) if ($p<$len-1);
    return $entropy - $e1 + $e2;
  }
  die;
}

sub change_op {
  my($str, $op) = @_;

  if ($op->[0] eq 'c') {
    my($p, $newc) = ($op->[1], $op->[2]);
    substr($$str, $p, 1) = $newc;
  } elsif ($op->[0] eq 'd') {
    my($p) = ($op->[1]);
    substr($$str, $p, 1) = "";
  } elsif ($op->[0] eq 'i') {
    my($p, $newc) = ($op->[1], $op->[2]);
    substr($$str, $p, 0) = $newc;
  } else {
    die;
  }
}

sub mutate {
  my($str, $num) = @_;

  my $entropy = entropy($str);
  while ($num>0) {
    my $len = length($str);

    my $op;

    my($Pdel) = (1-$Pchange)/(5+4/$len);
    my($Pins) = (1-$Pchange)*(4*$len+4)/(5*$len+4);

    die "Bad probs: $Pdel $Pins $Pchange" if (abs(1-$Pdel-$Pins-$Pchange)>0.0001);

    my $p = rand();
    if ($p<$Pchange) {
      # Try a change
      my $pos = int(rand($len));
      my $c = substr($str, $pos, 1);
      my $newc;
      do {
	$newc = (keys %$probs)[int(rand(keys %$probs))];
      } while ($newc eq $c);
      $op = ['c', $pos, $newc];
    } elsif ($p<$Pdel+$Pchange) {
      # Try a delete
      my $pos = int(rand($len));
      $op = ['d', $pos];
    } elsif ($p<$Pins+$Pdel+$Pchange) {
      # Try an insert
      my $pos = int(rand($len+1));
      my $newc = (keys %$probs)[int(rand(keys %$probs))];
      $op = ['i', $pos, $newc];
    } else {
      die;
    }


    my $entropy2 = mutate_entropy($entropy, \$str, $len, $op);

    my $r_ij = 2**($entropy-$entropy2);

    if ($r_ij > 1  ||		# Make change for certain
	rand() < $r_ij) {	# Make change probabilistically

      change_op(\$str, $op);
      $mutate_counts{$op->[0]}++;
      $entropy = $entropy2;

      $num--;
    }
  }
  return $str;
}




sub fit_gauss {
  my($data) = @_;
  my($sum,$sum_sq) = (0,0);
  for my $d (@$data) {
    $sum += $d;
    $sum_sq += $d*$d;
  }

  my $mean = $sum/@$data;
  my $var  = ($sum_sq/@$data) - ($mean*$mean);

  return ($mean, sqrt($var));
}

sub prob {
  my($hist, $c) = @_;
  return $probs->{substr($hist,-1)}{$c};
}

sub gen_sequence {		# Generate a sequence from 1st order markov
  my($len) = @_;

  my $str = (keys %$probs)[int(rand(keys %$probs))]; # Pick any old starting char
  for (2..$len) {
    my $last = substr($str, -1);
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

sub entropy {
  my($str) = @_;

  my $s = 0;
  for my $i (1..length($str)-1) {
    my $p += prob(substr($str,0,$i), substr($str,$i,1));
    $s += -log2($p);
  }
  return $s/length($str);
}

sub model_entropy {
  my $e = 0;
  for my $k1 (keys %$probs) {
    my $e2 = 0;
    for my $k2 (keys %$probs) {
      $e2 += -$probs->{$k1}{$k2}*log2($probs->{$k1}{$k2});
    }
    $e += $e2/4;
  }
  return $e;
}

sub model_power {
  my($pwr) = @_;
  for my $k1 (keys %$probs) {
    my $s;
    for my $k2 (keys %$probs) {
      $probs->{$k1}{$k2} **= $pwr;
      $s += $probs->{$k1}{$k2};
    }
    for my $k2 (keys %$probs) {
      $probs->{$k1}{$k2} /= $s;
    }
  }
}

sub rand_gauss {
  my($mu, $sd) = @_;

  # This method of generating gaussian distributed numbers is pretty bad, but it'll do.
  my $s = 0;
  for (1..12) { $s += rand(); }
  $s -= 6;
  $s *= $sd;
  $s += $mu;
  return $s;
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
