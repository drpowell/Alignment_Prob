
package Markov;

# Implement an arbitrary order adaptive Markov Model

use strict;

sub new {
  my($package, $order, $alphabet, $fade) = @_;
  my(@alphabet) = @{$alphabet};
  my($self) = {
    ORDER     => $order,
    ALPHABET  => [@alphabet],
    ALPHASIZE => scalar @alphabet,
    EVENTS => "",
    FADE      => $fade,
#    MML_ESTIMATE => 1,        # Use mml estimate for markov model
                              # Basically is (a+1/n)/(b+1)  (n=alpha size)
                              # instead of (a+1)/(b+n)
  };

  my(%counts);
#  %counts = map { $_ => { TOTAL=>0, COUNTS => {map {$_ => 0} @alphabet} } } 
#             genSub($order, @alphabet);
  # Now initialise counts
  $self->{COUNTS} = \%counts;

  return bless($self);
}

sub genSub {
  my($i, @alpha) = @_;

  ($i<=0) && (return "");

  my (@res);
  for $a (@alpha) {
    my(@r) = genSub($i-1, @alpha);
    push(@res, map { "$a$_" } @r)
  }
  return @res;
}


# Adds event 'a', updates counts, and returns probability of that event
# occuring
sub addEvent {
  my($self, $a) = @_;

  die "'$a' not in [@{$self->{ALPHABET}}]" 
      if (!grep {$a eq $_} @{$self->{ALPHABET}});

  my($order) = $self->{ORDER};
  my($len) = length($self->{EVENTS});

  $self->{EVENTS} .= "$a";
  return (1/$self->{ALPHASIZE}) if ($len < $order);

  my($past) = substr($self->{EVENTS}, $len-$order, $order);

  my($num) = $self->{COUNTS}{$past}{COUNTS}{$a} || 0;
  my($den) = $self->{COUNTS}{$past}{TOTAL} || 0;
  $self->{COUNTS}{$past}{TOTAL}++;
  $self->{COUNTS}{$past}{COUNTS}{$a}++;

  # Use total is = fade then halve all counts in that context
  if ($self->{FADE} && $self->{COUNTS}{$past}{TOTAL}==$self->{FADE}) {
    $self->fadeCounts($past);
  }

  if ($self->{MML_ESTIMATE}) {
    return ($num+(1/$self->{ALPHASIZE}))/($den+1);
  } else {
    return ($num+1)/($den+$self->{ALPHASIZE});
  }
}

sub fadeCounts {
  my($self, $past) = @_;
  my($tot) = 0;
  my($counts) = $self->{COUNTS}{$past}{COUNTS};
  for (keys %{$counts}) {
    $counts->{$_} = int($counts->{$_}/2);
    $counts->{$_} = 1 if (!$counts->{$_});
    $tot += $counts->{$_};
  }
  $self->{COUNTS}{$past}{TOTAL} = $tot;
}

# probEvent - probability of event 'a' occuring next
sub probEvent {
  my($self, $a, $past) = @_;

  die "'$a' not in [@{$self->{ALPHABET}}]" 
      if (!grep {$a eq $_} @{$self->{ALPHABET}});

  my($order) = $self->{ORDER};
  die "$past not of length $order" if (length($past) != $order);
  die "$past not in model ($past)" if (!defined($self->{COUNTS}{$past}{TOTAL}));
  
  my($num) = $self->{COUNTS}{$past}{COUNTS}{$a};
  my($den) = $self->{COUNTS}{$past}{TOTAL};

  if ($self->{MML_ESTIMATE}) {
    return ($num+(1/$self->{ALPHASIZE}))/($den+1);
  } else {
    return ($num+1)/($den+$self->{ALPHASIZE});
  }
}

# Check counts add up to their respective totals
sub sanityCheck {
  my($self) = @_;
  my($k,$a);
  for $k (keys %{$self->{COUNTS}}) {
    my($sum) = (0);
    $a = $self->{COUNTS}{$k};
    for (sort keys %{$a->{COUNTS}}) {
      $sum += $a->{COUNTS}{$_};
    }
    print "ERROR in context [$k] <$sum> != $a->{TOTAL}\n" if ($sum!=$a->{TOTAL});
  }
}

# displayCounts - display the MM counts
sub displayCounts {
  my($self,$fh) = @_;
#  $fh ||= fileno(STDOUT);

  my($k,$a);
  for $k (sort keys %{$self->{COUNTS}}) {
    $a = $self->{COUNTS}{$k};
    my($tot) = $a->{TOTAL};
    print "$k - $tot\n";
    for (sort keys %{$a->{COUNTS}}) {
      printf "  $_ = %6d  (%.4f)\n",
             $a->{COUNTS}{$_},$a->{COUNTS}{$_}/$tot;
    }
  }
}

# Generate all possible contexts for this order model
sub all_contexts {
  my($s, $o) = @_;
  (defined($o)) || ($o = $s->{ORDER});
  return [''] if ($o <=0 );

  my $t = $s->all_contexts($o-1);
  my $res = [];
  for my $a (@{$s->{ALPHABET}}) {
    for my $c (@$t) {
      push(@$res, $c . $a);
    }
  }
  return $res;
}


# model_entropy - Calculate the entropy of an arbitray order markov model.
# Note: Any high markov model can be considered as a first order
#       model with a large number of states
# Entropy = Sum_(over all contexts) ( P(context) * 
#                      Sum_(over all alpha chars) (P(a|context) * -log P(a|context) )
# To calculate this we need the steady-state state probabilities P(context).
# The transition probabilities can be used to form a number of linear equations:
#    P(context1)  = Sum_(over all contexts) ( P(context1| context) )
# With the equation Sum_(over all contexts) ( P(context) ) = 1
# These linear equations can be solved using Gaussian elimination.
sub model_entropy {
  my($s) = @_;
  return log2(scalar @{$s->{ALPHABET}}) if $s->{ORDER}<0;

  if ($s->{ORDER} >= 5) {
    print STDERR "Don't use model_entropy for models over order 4.  Too slow!\n";
    return log2(scalar @{$s->{ALPHABET}});
  }

  my %contexts;
  { my $i = 0;  map { $contexts{$_} = $i++; } @{ $s->all_contexts() }; }

#  print map { $_ . " -> " . $contexts{$_} . "\n" } keys %contexts;

  my($a_s) = ([]);
  for my $c (keys %contexts) {
    for my $k (@{$s->{ALPHABET}}) {
      my $to_c = ($s->{ORDER}==0 ? $k : substr($c, 1) . $k);
      my $to_i = ($s->{ORDER}==0 ? 0  : $contexts{$to_c});

      my($num,$den) = (0,0);
      if (exists($s->{COUNTS}{$c})) {
	$den = $s->{COUNTS}{$c}{TOTAL};
	$num = $s->{COUNTS}{$c}{COUNTS}{$k} || 0;
      }

      $a_s->[$to_i][$contexts{$c}] = ($num+1)/($den+$s->{ALPHASIZE});
    }
  }

  for my $i (0..scalar(@$a_s)-1) { $a_s->[$i][$i] -= 1; };

  my($b_s) = ([(0) x scalar(@$a_s)]);

  # Fill in all the undefined values with 0
  for my $r (@$a_s) {
    for my $i (0 .. scalar(@$a_s)-1) {
      $r->[$i] ||= 0;
    }
  }

  # Equations are under-constrained.  We can delete any 1 and replace it with sum(probs) = 1
  for my $j (0..scalar(@$a_s)-1) {
    $a_s->[0][$j] = 1;
  };
  $b_s->[0] = 1;

  solve_matrix($a_s, $b_s);

  # b_s now contains the equlibrium probabilities for each context
  my $entropy = 0;
  for my $c (keys %contexts) {
    my $h = 0;
    for my $k (@{$s->{ALPHABET}}) {
      my($num,$den) = (0,0);
      if (exists($s->{COUNTS}{$c})) {
	$den = $s->{COUNTS}{$c}{TOTAL};
	$num = $s->{COUNTS}{$c}{COUNTS}{$k} || 0;
      }

      my $p = ($num+1)/($den+$s->{ALPHASIZE});
      $h += - $p * log2($p);
    }
    $entropy += $h * $b_s->[$contexts{$c}];
  }

  return $entropy;
}


# ------------------------------------------------------------------------------------------

# Solve a set of linear equations of the form [As] * [xs] = [Bs]
sub solve_matrix {
  my($a_s, $b_s) = @_;

  my $DEBUG = 0;

  print_matrix($a_s,$b_s) if ($DEBUG);

  my($m,$n) = (scalar @$a_s, scalar @{$a_s->[0]});

  print "Doing $m rows, $n cols\n" if ($DEBUG);
  die "Only square matrices at the moment!" if ($m != $n);

  die "Bad rank in b_s" if $m != @$b_s;

  for my $i (0 .. $m-2) {
    # Find the row with the largest value in column[i]  (this will be the pivot)
    my $max = $i;
    for my $i2 ($i+1 .. $m-1) {
      $max = $i2 if (abs($a_s->[$i2][$i] > abs($a_s->[$max][$i])))
    }
    # Swap row $max with row $i  (and swap in the b vector)
    ($a_s->[$i], $a_s->[$max]) = ($a_s->[$max], $a_s->[$i]);
    ($b_s->[$i], $b_s->[$max]) = ($b_s->[$max], $b_s->[$i]);

    # row a_s[i] is what we are subtracting from everything

    for my $i2 ($i+1 .. $m-1) {
      # row a_s[i2] is the row we are about to subtract from
      next if ($a_s->[$i2][$i] == 0);
      my $r = $a_s->[$i2][$i] / $a_s->[$i][$i];

      for my $j (0 .. $n-1) {
	# j is the column
	$a_s->[$i2][$j] -= $a_s->[$i][$j] * $r;
      }
      $b_s->[$i2] -= $b_s->[$i] * $r;
    }
  }
  print_matrix($a_s,$b_s) if ($DEBUG);

  for (my $i=$m-1; $i>=0; $i--) {
    my $sum = 0;
    for my $j ($i+1 .. $m-1) {
      $sum += $b_s->[$j] * $a_s->[$i][$j];
      $a_s->[$i][$j] = 0;
    }

    if ($a_s->[$i][$i] == 0) {
      if ($b_s->[$i] == 0) {
	print "Equations are underconstrained.\n";
	last;
      } else {
	print "Equations are inconsistent.\n";
	last;
      }
    }

    $b_s->[$i] = ($b_s->[$i] - $sum) / $a_s->[$i][$i];
    $a_s->[$i][$i] = 1;
  }
  print_matrix($a_s,$b_s) if ($DEBUG);
}

sub print_matrix {
  my($a_s, $b_s) = @_;

  my($m,$n) = (scalar @$a_s, scalar @{$a_s->[0]});

  for my $i (0 .. $m - 1) {
    print "| ";
    for my $j (0 .. $n - 1) {
      printf "%5.2f ", $a_s->[$i][$j];
    }
    printf "| %s ", ($i == int($m/2) ? "=" : " ");
    printf " | %5.2f |\n", $b_s->[$i];
  }
  print "\n";

}

sub log2 {
  warn "Can't take log of 0" if ($_[0]==0);
  return log($_[0])/log(2); 
}

1;
