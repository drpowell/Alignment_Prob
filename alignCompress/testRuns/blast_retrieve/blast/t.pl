#!/usr/bin/perl -w

use Data::Dumper;

use Bio::Seq;
use Bio::Seq::RichSeq;


#my (@pos) = (571,872,969,1099,1219,1853);
my (@pos) = (571,1785);
#my (@pos) = (571,2085);
my $dat;

{ my $VAR1;
  my $str=join '', <>;
  $dat=eval($str);
}

for my $i (1..$#pos) {
  my $ss = $dat->trunc($pos[$i-1], $pos[$i]-1);
  printf "%4d - %4d :\n", $pos[$i-1], $pos[$i];
  my $str = $ss->revcom()->seq();
  my $pos = $pos[$i-1];
  while ($str =~ /(.{1,50})/g) {
    my $s = $1;
    $s =~ s/(.{1,10})/$1 /g;
    printf "    %s\n", $s;
#    printf "%4d  %s\n", $pos, $s;
    $pos = $pos[$i-1] + pos($str)
  }
  printf "\nTranslated: %s\n\n", $ss->revcom()->translate()->seq();
}

#print Dumper($ss),"\n";

