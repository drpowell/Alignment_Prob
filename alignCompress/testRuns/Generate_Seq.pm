
package Generate_Seq;

use strict;
use Data::Dumper;

my @alphabet = qw(a t g c);

sub new {
  my($p,$model,$params) = @_;

  my $mfunc;
  for ($model) {
    (/uniform/i) && do {
      *init  = \&uniform_init;
      *model = \&uniform_model;
      last;
    };

    (/markov0/i) && do {
      *init  = \&markov0_init;
      *model = \&markov0_model;
      last;
    };

    die "Unknown model '$model'";
  }

  my $s = {
	   MODEL   => $model,
	  };
  bless($s, $p);

  $s->init($params);

  return $s;
}

sub uniform_init {
  my($self, $alphabet) = @_;
  $self->{PARAMS} = $alphabet;
}
sub uniform_model {
  my($self) = @_;
  my $alphabet = $self->{PARAMS};
  return $alphabet[int(rand(@alphabet))];
}


sub markov0_init {
  my($self, $probs) = @_;
  my $c=0;
  for my $v (values %$probs) { $c += $v; }
  die "Probs sum to $c" if ($c != 1);
  $self->{PARAMS} = $probs;
}

sub markov0_model {
  my($self) = @_;
  my $probs = $self->{PARAMS};
  my $p = rand();
  for my $k (keys %$probs) {
    return $k if ($p < $probs->{$k});
    $p -= $probs->{$k};
  }
}


sub as_string {
  my($self) = @_;

  my $s;
  $s .= "Model=$self->{MODEL}\n";
  my $d = Data::Dumper->new([ $self->{PARAMS} ]);
  $s .= sprintf "  params=%s\n", $d->Terse(1)->Indent(2)->Dump();
}

sub unrelated {
  my($self, $len) = @_;

  my(@s1,@s2);
  for my $i (1..$len) {
    push(@s1, $self->model(\@s1));
    push(@s2, $self->model(\@s2));
  }
  return ((join '',@s1), join( '',@s2));
}

sub related {
  my($self, $len, $mutates) = @_;
  my(@s1,@s2);
  for my $i (1..$len) {
    push(@s1, $self->model(\@s1));
  }
  @s2 = @s1;

  for my $i (1..$mutates) {
    my $p = int(rand($#s2+1));
    $s2[$p] = $self->model([@s1[0..$p-1]]);
  }

  return ((join '',@s1), join( '',@s2));
}

1;

__DATA__
package main;

my $model = new Generate_Seq('markov0', {'a'=>0.1, 't'=>0.1, 'g'=>0.3, 'c'=>0.5});
#my $model = new Generate_Seq('uniform',[qw(a t g c)]);

#print join("\n",$model->unrelated(100)),"\n";
print join("\n",$model->related(40,2)),"\n";
