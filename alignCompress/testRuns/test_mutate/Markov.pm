
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

1;
