
package DPA;

use strict;

sub new {
  my($p, $s1, $s2, $match, $mis, $ins) = @_;

  if (!defined($match)) {
    $match = 0; $mis = 1; $ins = 1;
    #$match=2.2; $mis=7.81; $ins=6.71;   # These 2 set of costs should give the same rank order of alignments
    #$match=6; $mis=6.34; $ins=3.34;
  }

  my $s = {
	   S1 => $s1,
	   S2 => $s2,

	   MATCH => $match,
	   MIS   => $mis,
	   INDEL => $ins,
	  };
  bless($s, $p);

  $s->doDPA();

  return $s;
}

sub mis {
  my($self, $a,$b) = @_;
  return ($a eq $b ? $self->{MATCH} : $self->{MIS});
};

sub edit_dist {
  my($self) = @_;
  return $self->{EDIT_DIST};
}

sub doDPA {
  my($self) = @_;
  my @a = split(//,$self->{S1});
  my @b = split(//,$self->{S2});

  my($D,$P) = ([],[]);
  for my $i (0..@a) {
    for my $j (0..@b) {
      ($i==0 && $j==0) && do {$D->[0][0] = 0; next};
      ($i==0) && do { $D->[0][$j] = $D->[0][$j-1] + $self->{INDEL}; $P->[0][$j]=2; next};
      ($j==0) && do { $D->[$i][0] = $D->[$i-1][0] + $self->{INDEL}; $P->[$i][0]=1; next};

      $D->[$i][$j] = min( $D->[$i-1][$j-1] + $self->mis($a[$i-1],$b[$j-1]),
			$D->[$i-1][$j] + $self->{INDEL},
			$D->[$i][$j-1] + $self->{INDEL});

      $P->[$i][$j] = min_( $D->[$i-1][$j-1] + $self->mis($a[$i-1],$b[$j-1]),
			 $D->[$i-1][$j] + $self->{INDEL},
			 $D->[$i][$j-1] + $self->{INDEL});
    }
  }

  $self->{D} = $D;
  $self->{P} = $P;
  $self->{EDIT_DIST} = $D->[-1][-1];
}

sub getAlignment {
  my($self) = @_;
  my @a = split(//,$self->{S1});
  my @b = split(//,$self->{S2});
  my($i,$j) = (scalar @a, scalar @b);
  my($ra,$rb) = ("","");
  while ($i!=0 || $j!=0) {
#    $D->[$i][$j]='*';     # Put '*' on the optimal alignment
    ($self->{P}[$i][$j] == 0) && do {$i--; $j--; $ra .= $a[$i]; $rb .= $b[$j]; next};
    ($self->{P}[$i][$j] == 1) && do {$i--; $ra.= $a[$i]; $rb .= '-'; next};
    ($self->{P}[$i][$j] == 2) && do {$j--; $ra .= '-'; $rb .= $b[$j];next};
  }
  $ra = join('',reverse split(//,$ra));
  $rb = join('',reverse split(//,$rb));
  for my $i (0..length($ra)-1) {
    my $s = substr($ra,$i,1) eq substr($rb,$i,1);
    substr($ra,$i,1) = ($s ? uc(substr($ra,$i,1)) : lc(substr($ra,$i,1)));
    substr($rb,$i,1) = ($s ? uc(substr($rb,$i,1)) : lc(substr($rb,$i,1)));
  }

  return ($ra,$rb);
}

sub printMatrix {
  my($self) = @_;
  for my $i (0..length($self->{S1})) {
    for my $j (0..length($self->{S2})) {
      printf "%3s ",($self->{D}[$i][$j]);
    }
    printf "\n";
  }
}

sub min_ {
  my($a,$b,$c) = @_;
  
  ($a<$b) ? (($a<$c) ? 0 : 2) : (($b<$c) ? 1 : 2);
};

sub min {
  my($a,$b,$c) = @_;
  
  ($a<$b) ? (($a<$c) ? $a : $c) : (($b<$c) ? $b : $c);
};

1;

__DATA__

package main;

my $A=<>;
my $B=<>;

chomp($A);
chomp($B);

my $d = new DPA($A,$B);

print "Edit distance = " . $d->edit_dist() . "\n";


my($ra,$rb) = $d->getAlignment();
print "$ra\n$rb\n";
$d->printMatrix();

