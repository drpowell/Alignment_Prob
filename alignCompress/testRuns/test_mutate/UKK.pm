
package UKK;

use strict;

our $INDEL = 1;
our $MISMATCH = 1;

sub new {
  my($p, $s1, $s2) = @_;

  my $s = {
	   A => [split(//,$s1)],
	   B => [split(//,$s2)],
	   U => [],
	   COMPUTED => [],
	   D => [],
	   abOffset => length($s2),
	  };
  bless($s,$p);

  $s->doUkk();

  return $s;
}

sub doUkk {
  my($self) = @_;

  my $finalab = (scalar @{$self->{A}}) - (scalar @{$self->{B}});
  my $cost=-1;
  my $dist;
  do {
    $cost++;
    $dist = $self->Ukk($finalab,$cost);
  } while ($dist < @{$self->{A}});

  $self->{COST} = $cost;
}


sub edit_dist {
  my($self) = @_;
  return $self->{COST};
}

sub Ukk {
  my($self, $ab, $cost) = @_;

  return -1 if ($ab==0 && $cost==-1);

  return -1E10 if (abs($ab)>$cost);

  if (defined($self->{COMPUTED}[$ab+$self->{abOffset}][$cost]) && 
              $self->{COMPUTED}[$ab+$self->{abOffset}][$cost]) {
    return $self->{U}[$ab+$self->{abOffset}][$cost];
  }

  my($v1,$v2,$v3) = ($self->Ukk($ab+1, $cost-$INDEL),
                     $self->Ukk($ab,   $cost-$MISMATCH)+1,
                     $self->Ukk($ab-1, $cost-$INDEL)+1);

  my $dist = max($v1,$v2,$v3);

  while (defined($self->{A}[$dist]) && defined($self->{B}[$dist-$ab]) &&
         $self->{A}[$dist] eq $self->{B}[$dist-$ab]) {
    $self->{D}[$dist][$dist-$ab] = ":";
    $dist++;
  }

  $self->{U}[$ab+$self->{abOffset}][$cost] = $dist;
  $self->{COMPUTED}[$ab+$self->{abOffset}][$cost] = 1;

  $self->{D}[$dist][$dist-$ab] = $cost;

  return $dist;
}


sub max {
  my($i);
  for my $v (@_) {
    (!defined($i) || $v>$i) && ($i=$v);
  }
  return $i;
}

1;

__DATA__

package main;

my $A=<>;
my $B=<>;

chomp($A);
chomp($B);

my $d = new UKK($A,$B);

print "Edit distance = " . $d->edit_dist() . "\n";

__DATA__

sub printDPAmatrix {
  print "DPA matrix\n";
  for my $i (0..@A) {
    for my $j (0..@B) {
      printf "%3s ",(defined($D[$i][$j]) ? $D[$i][$j] : "-");
    }
    printf "\n";
  }
}

sub printUKKmatrix {
  my($cost) = @_;
  print "UKK matrix for cost=$cost\n";
  for my $ab (-$cost .. $cost) {
    printf "%3d: ",$ab;
    for my $d (0..$cost) {
#      if (defined($U[$ab + $abOffset][$d])) {
      if (defined($computed[$ab + $abOffset][$d]) &&
          $computed[$ab + $abOffset][$d]) {
        printf "%2d ",$U[$ab + $abOffset][$d];
      } else {
        print "   ";
      }
    }
    print "\n";
  }
}

