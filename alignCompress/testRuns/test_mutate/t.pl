#!/usr/bin/perl -w


# This script attempts to mutate a sequence generated from a 1st order
# Markov model while still keeping the mutated sequence _in_ the model.
# This approach is inspired by discussion with Lloyd and Wallace's proof
# for the Metropolis algorithm (http://www.csse.monash.edu.au/~lloyd/tildeMML/Local/1995-Metropolis/)
#
# Works as follows:
#    Generate sequence s1 using MM.
#    Determine number of possible mutations.
#         - each position has 3 possible changes (using DNA alphabet)
#         - each position has 1 possible delete
#         - between each position has 4 possible inserts
#    This is the number of out arcs, di, for state Si
#    From these mutations pick one uniformly randomly, it leads to state Sj
#    Calculate dj in a similar manner to di
#    Calculate ratio of Pi/Pj
#    Calculate R_ij = (Pj/Pi) * (di/dj)
#    if (R_ij > 1) then make the mutation to s1 that corresponds to moving to state Sj
#    if (R_ij < 1) then make the mutation with probability = R_ij
#
# Since we don't take sequence length into account, short sequence should be favoured.
# But on the other hand, there is one more insertation point than deletion, and every insertion
# position has 4 possible characters, so many more insertion arcs.
# These may have a tendency to balance each other out?




use strict;
use Data::Dumper;

use UKK;

use vars qw($main_loop);	# Global for main loop counter. So can access in signal handler
$SIG{HUP} = sub {print STDERR "Upto: $main_loop\n";};


#my($len_lower, $len_upper) = (100,200);
my($len_lower, $len_upper) = (100,100);

my $probs = { 'a' => {'a' => 0.11, 't'=>0.09, 'g'=>0.3,  'c'=>0.5},
	      't' => {'a' => 0.09, 't'=>0.11, 'g'=>0.5,  'c'=>0.3},
	      'g' => {'a' => 0.5,  't'=>0.3,  'g'=>0.09, 'c'=>0.11},
	      'c' => {'a' => 0.5,  't'=>0.11, 'g'=>0.09, 'c'=>0.3}
	    };

model_power(0.001);
printf "Model = %s\n", Dumper($probs);
printf "Model entropy = %f\n", model_entropy();

my(%mutate_counts);

my(@e1,@e2);
my(@l1,@l2);
my(@ed);
for $main_loop (1..1000) {
  my $str = gen_sequence(int(rand($len_upper-$len_lower+1))+$len_lower);
  push(@e1, entropy($str));
  push(@l1, length($str));

#  my $str2 = gen_sequence(int(rand($len_upper-$len_lower+1))+$len_lower);
  my $str2 = mutate($str, 1);
  push(@e2, entropy($str2));
  push(@l2, length($str2));

#  my $ukk = new UKK($str,$str2);
#  push(@ed, $ukk->edit_dist);

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

sub possible_changes {
  my($len) = @_;
  return $len*3;
}

sub possible_deletes {
  my($len) = @_;
  return $len;
}

sub possible_inserts {
  my($len) = @_;
  return ($len+1)*4;
}

sub possible_mutates {		# This assumes DNA alphabet
  my($len)=@_;
  return (possible_changes($len)+
	  possible_deletes($len)+
	  possible_inserts($len));
}

sub possible_mutates2 {		# This assumes DNA alphabet
  my($op, $len)=@_;
  ($op->[0] eq 'i') && ($len++);
  ($op->[0] eq 'd') && ($len--);
  return possible_mutates($len);
}

sub uniform_mutate {
  my($str,$len) = @_;
  my $n = possible_mutates($len);
  my $p = int(rand($n));
  if ($p < possible_changes($len)) {
    my $c = $p%3;
    $p = int($p/3);

    my @chars;
    for my $k (keys %$probs) { push @chars, $k if $k ne substr($$str, $p, 1); };
    return ('c', $p, $chars[$c]);
  } else {
    $p -= possible_changes($len);
    if ($p < possible_deletes($len)) {
      return ('d', $p);
    } else {
      $p -= possible_deletes($len);
      ($p < possible_inserts($len)) || die "Bad rand num $p";
      my $c = $p%4;
      $p = int($p/4);
      return ('i', $p, (keys %$probs)[$c]);
    }
  }
  die;
}

sub mutate_entropy {
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
    my $di = possible_mutates($len);

    my @op = uniform_mutate(\$str,$len);

    my $dj = possible_mutates2(\@op,$len);
    my $entropy2 = mutate_entropy($entropy, \$str, $len, \@op);

    my $r_ij = $di/$dj * (2**($entropy-$entropy2));

    if ($r_ij > 1  ||		# Make change for certain
	rand() < $r_ij) {	# Make change probabilistically

      change_op(\$str, \@op);
      $mutate_counts{$op[0]}++;
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
      $p -= prob($str, $cc);
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

  my $s = log2(scalar keys %$probs); # Uniform dist on first character
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
