
package Multi_model;

use strict;
use Data::Dumper;
use Markov_gen;


{
  my($me) = ($0 =~ m{([^/]*)\.pm$});
  if ($me && $me eq __PACKAGE__) {
    $^W = 1;
    MAIN_TEST(@ARGV);
  }
}

sub MAIN_TEST {
  my $m = new Multi_model();
  my($seq, $entropy) = ($m->gen_sequence(200));
  print "$seq\n";
  print "Entropy (if cut-points known) = ",$entropy," bits/char\n";
  $seq = lc($seq);

  my $len = $m->fit_model($seq);
  printf "Encoding len = $len    entropy (per char)=%f\n", $len/length($seq);
}



sub new {
  my $p = shift;
  my $self  = {
	       SWITCH_PROB => 0.5,
	       SUB_LEN     => 20,
	       SUB_RANGE   => 10,

	       M1 => new Markov_gen(-1, [qw(a t g c)]),

	       M2 => new Markov_gen(1, [qw(a t g c)],
				    {
				     'a' => {a=>0.28, t=>0.7,  g=>0.01, c=>0.01},
				     't' => {a=>0.7,  t=>0.28, g=>0.01, c=>0.01},
				     'g' => {a=>0.49, t=>0.49, g=>0.01, c=>0.01},
				     'c' => {a=>0.49, t=>0.49, g=>0.01, c=>0.01}
				    }),
	      };

  bless($self, $p);
  return $self;
}

sub as_string {
  my($self) = @_;

  my $s = "Blend Model: ";
  $s .= sprintf("Pswitch=%f SUB_LEN=%d SUB_RANGE=%d\n", 
		$self->{SWITCH_PROB}, $self->{SUB_LEN}, $self->{SUB_RANGE});
  $s .= "    Model1: " . $self->{M1}->as_string();
  $s .= "    Model2: " . $self->{M2}->as_string();
  return $s;
}

sub gen_sequence {
  my($self, $len) = @_;

  my $entropy=0;
  my $seq;
  while (1) {
    my $m = rand()<$self->{SWITCH_PROB} ? $self->{M1} : $self->{M2};
    my $s = $m->gen_sequence(rand_length($self->{SUB_LEN}, $self->{SUB_RANGE}));

    $entropy += $m->entropy($s)*length($s);

    ($m == $self->{M2}) && ($s=uc($s));

    $seq .= $s;

    last if length($seq)>=$len;
  }

  if (wantarray) {
    return ($seq, $entropy/length($seq));
  } else {
    return $seq;
  }
}

# Just use one of the models to do the mutations.
# This is far from ideal, but doing it properly would be difficult
sub mutate {
  my($self, $seq, $numMutations) = @_;
  return $self->{M1}->mutate($seq, $numMutations);
}

sub predict {
  my($model, $seq) = @_;
  my $probs = $model->{PROBS};
  my $alpha = $model->{ALPHA};
  my $order = ($model->{ORDER} < 0 ? 0 : $model->{ORDER});

  if (length($seq) < $order) {
    return { map { $_ => 1/(scalar @$alpha) } @$alpha };
  }

  my $last = substr($seq, length($seq)-$order);
  return { %{$probs->{$last}} };
}

sub fit_model {
  my($self, $seq) = @_;
  my(@models) = ($self->{M1}, $self->{M2});
  my(@m_lens) = (1) x @models;

  my $len = 0;
  for my $i (0 .. length($seq)-1) {
    my @weights;
    my $sum = 1E20;
    map { $sum = logplus($sum, $_) } @m_lens;
    for my $mi (0..$#models) {
      $weights[$mi] = exp2($sum - $m_lens[$mi]);
    }

    my $p = 0;
    my $char = substr($seq, $i, 1);
    for my $mi (0..$#models) {
      my $pred = predict($models[$mi], substr($seq, 0, $i));

      $p += $pred->{$char} * $weights[$mi];

      $m_lens[$mi] += -log2( $pred->{$char} );

      # Fade
      $m_lens[$mi] /= 1.8;
    }

    $len += -log2( $p );
  }

  return $len;
}






sub rand_length {
  my($mid, $range) = @_;
  return $mid + int(rand($range*2+1))-$range;
}

sub log2 {
  warn "Can't take log of 0" if ($_[0]==0);
  return log($_[0])/log(2); 
}

sub exp2 {
  return exp($_[0] * log(2));
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

sub logplus {
  my($a,$b) = @_;
  if ($a>$b) {
    my $t = $a;
    $a = $b;
    $b = $t;
  }
  if ($b>$a+30) { return $a }
  return $a - log2(1+exp2($a-$b));
}


1;
