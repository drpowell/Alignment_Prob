
package Markov_gen;

# Can be used to generate 0, 1st or 2nd order Markov models and sequence from the models.
# Use randModel() to generate a random model.
# If a more uniform model is desired use model_power($pwr) with $pwr<1  (use normalise() afterwards)
# If a more extreme model is desired use model_power($pwr) with $pwr>1  (use normalise() afterwards)

# You can also make a 1st order model that has uniform 0 order statistics: makeModel($pwr)
# You can also make a 2nd order model that has uniform 1st order statistics: makeModel($pwr)
# Orginally this was to show false positives in the fasta S-W algorithm

# The method negate_model() changes each probability, p, to 1-p.  Then normalises.
# This can be used to generate sequence that have high entropy under some model.

# The method copy_model() returns a deep copy of the model object.

# The method mutate($str,$num) will make $num mutations to $str.  This will attempt to
# keep the result matching the Markov model.  This is done by making mutation probabilistically
# based on how the mutation changes the entropy of the sequence.
# The probability $self->{PCHANGE} may be set to influence the ratio of changes to indels.

use strict;
use Data::Dumper;

{
  my($me) = ($0 =~ m{([^/]*)\.pm$});
  if ($me && $me eq __PACKAGE__) {
    $^W = 1;
    MAIN_TEST(@ARGV);
  }
}


sub MAIN_TEST {
  my $order = 1;
  my $t = new Markov_gen($order, [qw(a t g c)], 
			 { 'a' => { 'a' => 11, 't' => 2, 'g'=> 1, 'c' => 1},
			   't' => { 'a' => 1,  't' => 12,'g'=> 1, 'c' => 1},
			   'g' => { 'a' => 1,  't' => 1, 'g'=> 11,'c' => 2},
			   'c' => { 'a' => 1,  't' => 1, 'g'=> 1, 'c' => 12}
			 });
  $t->normalise();

  print $t->as_string(),"\n";
  print "New Entropy = ",$t->model_entropy(),"\n";
  print "Old Entropy = ",$t->model_entropy_wrong(),"\n";

  my $tot=0;
  for (1..5) {
    my $seq = $t->gen_sequence(100000);
    my $v = $t->entropy($seq);
    $tot += $v;
    print "A sequence entropy = ",$v,"\n";
  }
  print "Avg = ", $tot/5,"\n";

  #my $s1 = $t->gen_sequence(100);
  #my $s2 = $t->mutate($s1, 50);
  #print $t->gen_sequence(10) . $s1 . "\n";
  #print $s2 . $t->gen_sequence(10) . "\n";
  #my $t = new Markov_gen(0, [qw(a t g c)]);
  #$t->model_power(3);
  #$t->makeUniModel(5);
  #my $t2 = $t->copy_model();
  #$t2->negate_model();

  #print $t->as_string(),"\n";
  #print $t2->as_string(),"\n";

  #my $s1 = $t->gen_sequence(10);
  #my $s2 = $t2->gen_sequence(10);

  #print "s1=$s1\ns2=$s2\n";

  #printf "s1 under t  = %f\n",$t->entropy($s1);
  #printf "s1 under t2 = %f\n",$t2->entropy($s1);
  #printf "s2 under t  = %f\n",$t->entropy($s2);
  #printf "s2 under t2 = %f\n",$t2->entropy($s2);
}

sub new {
  my($p, $order, $alphabet, $probs) = @_;

  if ($order<-1 || $order>2) {
    die "Cannot do order $order models";
  }

  if (ref($alphabet) ne 'ARRAY') {
    die "<alphabet> must be an array ref";
  }


  my $s = {
	   ORDER => $order,
	   ALPHA => [@$alphabet],
	   PROBS => $probs,

	   PCHANGE => 0.6,	# Used by the mutate() func

	  };
  bless($s,$p);

  if (!defined($s->{PROBS})) { $s->randModel(); }

  return $s;
}

sub randModel {
  my($s) = @_;
  ($s->{ORDER} == -1) && (return randModel_order_uniform(@_));
  ($s->{ORDER} ==  0) && (return randModel_order0(@_));
  ($s->{ORDER} ==  1) && (return randModel_order1(@_));
  ($s->{ORDER} ==  2) && (return randModel_order2(@_));
  die "Bad order";
}

sub fixStats {
  my($s) = @_;
  ($s->{ORDER} == -1) && (die "fixStats - not implemented for -1th order");
  ($s->{ORDER} ==  0) && (die "fixStats - not implemented for 0th order");
  ($s->{ORDER} ==  1) && (return fixStats_order1(@_));
  ($s->{ORDER} ==  2) && (return fixStats_order2(@_));
  die "Bad order";
}

sub mutate_entropy {
  my($s) = @_;
  ($s->{ORDER} == -1) && (return mutate_entropy_order0(@_));
  ($s->{ORDER} ==  0) && (return mutate_entropy_order0(@_));
  ($s->{ORDER} ==  1) && (return mutate_entropy_order1(@_));
  ($s->{ORDER} ==  2) && do {
    print STDERR "mutate_entropy not implemented for 2nd order. Using 0-order mutate\n";
    return mutate_entropy_order_uni(@_);
  };
    
  die "Bad order";
}

sub as_string {
  my($self) = @_;

  my $s;
  $s .= "Model order=$self->{ORDER}";
  $s .= sprintf " entropy=%.4f", $self->model_entropy();
  $s .= sprintf " Pchange=%.2f", $self->{PCHANGE};
  my $d = Data::Dumper->new([ $self->{PROBS} ]);
  $s .= sprintf "  probs=%s\n", $d->Terse(1)->Indent(2)->Dump();
}

sub makeUniModel {
  my($s, $power) = @_;

  (!defined($power)) && ($power=1);

  $s->randModel();
  while ($s->fixStats() > 1E-10) { };
  if ($power != 1) {
    $s->model_power($power);
    while ($s->fixStats() > 1E-10) { };
  }
}

# Not random at all, uniform!
sub randModel_order_uniform {
  my($s) = @_;
  my $m = {};
  for my $i (@{$s->{ALPHA}}) { $m->{''}{$i} = 1 };
  $s->{PROBS} = $m;
  $s->normalise();
  return $m;
}

sub randModel_order0 {
  my($s) = @_;
  my $m = {};
  for my $i (@{$s->{ALPHA}}) {
	$m->{''}{$i} = rand();
  }
  $s->{PROBS} = $m;
  $s->normalise();
  return $m;
}

sub randModel_order1 {
  my($s) = @_;
  my $m = {};
  for my $i (@{$s->{ALPHA}}) {
    for my $j (@{$s->{ALPHA}}) {
	$m->{$i}{$j} = rand();
    }
  }
  $s->{PROBS} = $m;
  $s->normalise();
  return $m;
}

sub randModel_order2 {
  my($s) = @_;
  my $m = {};
  for my $i (@{$s->{ALPHA}}) {
    for my $j (@{$s->{ALPHA}}) {
      for my $k (@{$s->{ALPHA}}) {
	$m->{$i.$j}{$k} = rand();
      }
    }
  }
  $s->{PROBS} = $m;
  $s->normalise();
  return $m;
}

sub fixStats_order1 {
  my($s) = @_;
  my $m = $s->{PROBS};
  my($max) = $s->normalise();

  for my $k (@{$s->{ALPHA}}) {
    my $sum = 0;
    for my $i (@{$s->{ALPHA}}) {	
      $sum += $m->{$i}{$k};
    }
    for my $i (@{$s->{ALPHA}}) {	
      $m->{$i}{$k} /= $sum;
    }
    $max = max($max, abs($sum-1));
  }

  return $max;
}

sub fixStats_order2 {
  my($s) = @_;
  my $m = $s->{PROBS};
  my($max) = $s->normalise();

  for my $j (@{$s->{ALPHA}}) {
    for my $k (@{$s->{ALPHA}}) {
      my $sum = 0;
      for my $i (@{$s->{ALPHA}}) {	
	$sum += $m->{$i.$j}{$k};
      }
      for my $i (@{$s->{ALPHA}}) {	
	$m->{$i.$j}{$k} /= $sum;
      }
      $max = max($max, abs($sum-1));
    }
  }

  return $max;
}

sub normalise {
  my($s) = @_;
  my $m = $s->{PROBS};
  my $max=0;
  for my $k1 (keys %$m) {
    my $sum = 0;
    for my $k2 (keys %{$m->{$k1}}) {
      $sum += $m->{$k1}{$k2};
    }
    for my $k2 (keys %{$m->{$k1}}) {
      $m->{$k1}{$k2} /= $sum;
    }
    $max = max($max, abs($sum-1));
  }
  return $max;
}

# Calculate the entropy of an arbitray order markov model.
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
  return log2(scalar @{$s->{ALPHA}}) if $s->{ORDER}<0;

  my $m = $s->{PROBS};
  my %contexts;
  { my $i = 0;  map { $contexts{$_} = $i++; } keys %$m; }
  #%contexts = (a => 0, t=>1, g=>2, c=>3);

  print map { $_ . " -> " . $contexts{$_} . "\n" } keys %contexts;

  my($a_s) = ([]);
  for my $c (keys %contexts) {
    for my $k (keys %{$m->{$c}}) {
      my $to_c = substr($c, 1) . $k;
      my $to_i = ($s->{ORDER}==0 ? 0 : $contexts{$to_c});
      $a_s->[$to_i][$contexts{$c}] = $m->{$c}{$k};
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
    for my $k (keys %{$m->{$c}}) {
      $h += - $m->{$c}{$k} * log2($m->{$c}{$k});
    }
    $entropy += $h * $b_s->[$contexts{$c}];
  }

  return $entropy;
}

# model_entropy - calculate the entropy of the model (I think this is right)
#     THIS IS WRONG. OFTEN CLOSE BUT WRONG
sub model_entropy_wrong {
  my($s) = @_;
  my $m = $s->{PROBS};
  my $sum=0;
  for my $k1 (keys %$m) {
    for my $k2 (keys %{$m->{$k1}}) {
      $sum += $m->{$k1}{$k2} * -log2($m->{$k1}{$k2});
    }
  }
  $sum /= (scalar keys %$m);
  return $sum;
}

sub model_power {
  my($s, $pwr) = @_;
  return if ($s->{ORDER}<0);
  my $probs = $s->{PROBS};
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

sub copy_model {
  my($m) = @_;
  my $r;
  if (ref($m) eq 'Markov_gen') {
    for my $k (keys %$m) { $r->{copy_model($k)} = copy_model($m->{$k}) };
    bless($r, ref($m));
  } elsif (ref($m) eq 'HASH') {
    for my $k (keys %$m) { $r->{copy_model($k)} = copy_model($m->{$k}) };
  } elsif (ref($m) eq 'ARRAY') {
    for my $i (0..$#$m) { $r->[$i] = copy_model($m->[$i]); };
  } elsif (ref($m) eq 'SCALAR') {
    $r = $$m;
    $r = \$r;
  } elsif (ref($m) eq '') {
    $r = $m;
  } else {
    die "Don't know how to copy a ".ref($m)."\n";
  }
  return $r;
}

sub negate_model {
  my($s) = @_;
  return if ($s->{ORDER}<0);
  my $probs = $s->{PROBS};
  for my $k1 (keys %$probs) {
    for my $k2 (keys %{$probs->{$k1}}) {
      $probs->{$k1}{$k2} = 1-$probs->{$k1}{$k2};
    }
  }
  $s->normalise();
}

# Calculate the entropy of $str under the current model
sub entropy {
  my($s, $str) = @_;

  my $probs = $s->{PROBS};

  return undef if (length($str)==0);
  $str = lc $str;

  my $order = ($s->{ORDER} < 0 ? 0 : $s->{ORDER});

  my $sum = 0;
  for my $i (1 .. $order) {
    $sum += log2(scalar @{$s->{ALPHA}}); # Uniform for first chars
  }

  for my $i ($order .. length($str)-1) {
    my $hist = substr($str, $i-$order, $order);
    my $p = $probs->{$hist}{substr($str,$i,1)};
    $sum += -log2($p);
  }
  return $sum/length($str);
}

sub gen_sequence {		# Generate a sequence 
  my($self, $len) = @_;

  my $probs = $self->{PROBS};
  my $alpha = $self->{ALPHA};
  my $order = ($self->{ORDER} < 0 ? 0 : $self->{ORDER});

  my $str = ""; # Pick any $order starting chars
  for (1 .. $order) {
    $str .= $alpha->[rand(@$alpha)];
  }

  for (length($str)+1 .. $len) {
    my $last = substr($str, length($str)-$order);
    my $p = rand();
    my $c;
    for my $cc (keys %{$probs->{$last}}) {
      $p -= $probs->{$last}{$cc};
      if ($p <= 0) {
	$c = $cc;
	last;
      }
    }
    die "No char decided upon p=$p" if !defined($c);
    $str .= $c;
  }
  return $str;
}


#----------------------------------------------------------------------------
# The mutate functions

# These functions attempt to mutate a sequence generated from a
# Markov model while still keeping the mutated sequence _in_ the model.
# This approach is inspired by discussion with Lloyd and Wallace's proof
# for the Metropolis algorithm (http://www.csse.monash.edu.au/~lloyd/tildeMML/Local/1995-Metropolis/)
#
# Works as follows:
#    Generate sequence s1 using MM.
#



# Entropy doesn't change under a uniform model (I think)
sub mutate_entropy_order_uni {
  my($s, $entropy, $str, $len, $op) = @_;
  return $entropy;
}

sub mutate_entropy_order0 {
  # This is for 0th order MMs only!
  my($s, $entropy, $str, $len, $op) = @_;

  my $probs = $s->{PROBS};

  if ($op->[0] eq 'c') {
    my($p, $newc) = ($op->[1], $op->[2]);
    my $cc = lc substr($$str, $p,   1);
    my $e1 = -log2($probs->{''}{$cc});
    my $e2 = -log2($probs->{''}{$newc});
    return $entropy - $e1 + $e2;
  }

  if ($op->[0] eq 'd') {
    my($p) = ($op->[1]);
    my $cc = lc substr($$str, $p,   1);
    my $e1 = -log2($probs->{''}{$cc});
    return $entropy - $e1;
  }

  if ($op->[0] eq 'i') {
    my($p, $newc) = ($op->[1], $op->[2]);
    my $e2 = -log2($probs->{''}{$newc});
    return $entropy + $e2;
  }

  die;
}

sub mutate_entropy_order1 {
  # This is for 1st order MMs only!
  my($s, $entropy, $str, $len, $op) = @_;

  my $probs = $s->{PROBS};

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

# mutate - mutate 'str', 'num' times, but try to keep the result _in_ the model
sub mutate {
  my($s, $str, $num) = @_;

  my $probs = $s->{PROBS};

  my $entropy = $s->entropy($str);
  my %mutate_counts;
  while ($num>0) {
    my $len = length($str);

    my $op;

    my($Pchange, $Pdel, $Pins);
    $Pchange = ($len==0 ? 0 : $s->{PCHANGE});
    $Pdel    = ($len==0 ? 0 : (1-$Pchange)/(5+4/$len));
    $Pins    = ($len==0 ? 1 : (1-$Pchange)*(4*$len+4)/(5*$len+4));

    die "Bad probs: $Pdel $Pins $Pchange" if (abs(1-$Pdel-$Pins-$Pchange)>0.0001);

    my $p = rand();
    if ($p<$Pchange) {
      # Try a change
      my $pos = int(rand($len));
      my $c = substr($str, $pos, 1);
      my $newc;
      do {
	$newc = $s->{ALPHA}[int(rand(@{$s->{ALPHA}}))];
      } while ($newc eq $c);
      $op = ['c', $pos, $newc];
    } elsif ($p<$Pdel+$Pchange) {
      # Try a delete
      my $pos = int(rand($len));
      $op = ['d', $pos];
    } elsif ($p<$Pins+$Pdel+$Pchange) {
      # Try an insert
      my $pos = int(rand($len+1));
      my $newc = $s->{ALPHA}[int(rand(@{$s->{ALPHA}}))];
      $op = ['i', $pos, $newc];
    } else {
      die;
    }


    my $entropy2 = $s->mutate_entropy($entropy, \$str, $len, $op);

    my $r_ij = 2**($entropy-$entropy2);

    if ($r_ij > 1  ||		# Make change for certain
	rand() < $r_ij) {	# Make change probabilistically

      change_op(\$str, $op);
      $mutate_counts{$op->[0]}++;
      $entropy = $entropy2;

      $num--;
    }
  }

  if (wantarray) {
    return ($str, \%mutate_counts);
  } else {
    return $str;
  }
}


#----------------------------------------------------------------------------
# Misc funcs

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

sub max {
  my($i) = pop;
  for (@_) { ($_>$i) && ($i=$_); }
  $i;
}

sub min {
  my($i) = pop;
  for (@_) { ($_<$i) && ($i=$_); }
  $i;
}


# Solve a set of linear equations of the form [As] * [xs] = [Bs]
sub solve_matrix {
  my($a_s, $b_s) = @_;

  my $DEBUG = 1;

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




1;
