#!/usr/bin/perl -w
#
# Calculate the entropy under Markov Models of different orders

use Getopt::Long;
use Markov;
#use ReadGenBankFmt;

$| = 1;     # Flush stdout every print

my($lowOrder,$highOrder,$fade,$dispCounts);

my @Options = (
  {OPT=>"help",      VAR=>\&usage,      DESC=>"This help"},
  {OPT=>"fade=i",    VAR=>\$fade, DESC=>"Fade value (0 - nofade)", DEFAULT=>0},
  {OPT=>"start=i",   VAR=>\$lowOrder,   DESC=>"Starting order MM", DEFAULT=>0},
  {OPT=>"end=i",     VAR=>\$highOrder,  DESC=>"Ending order MM",   DEFAULT=>5},
  {OPT=>"counts!",   VAR=>\$dispCounts, DESC=>"Display MM counts", DEFAULT=>0},

);

&GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

# Now setup default values.
foreach (@Options) {
  if (defined($_->{DEFAULT}) && !${$_->{VAR}}) {
    ${$_->{VAR}} = $_->{DEFAULT};
  }
}

$str = <>; chomp($str);
@str = split(//,uc($str));

my @alpha = qw(A T G C);

printf "Alphabet size = %d  Null bits per char=%.2f bits\n",
        scalar @alpha, log(@alpha)/log(2);
printf "Len=%d fade=%d\n", scalar @str,$fade;

$null = (scalar @str)*log(@alpha)/log(2);
printf "Null encoding = %.2f bits\n", $null;

for $order ($lowOrder..$highOrder) {
print STDERR "Initializing $order MM\n";
  $m = new Markov($order, [@alpha], $fade);
print STDERR "Starting $order MM\n";
  my($len) = 0;
  for (@str) {
    $p = $m->addEvent($_);
    $len += -log($p);
  }
  printf("($order-order MM) Len = %.2f bits (Saving = %.2f bits) (%.4f bits/char) (model = %.4f bits/char)\n",
	 $len/log(2), $null-$len/log(2), $len/log(2)/(scalar @str),
	 $m->model_entropy());

  $m->sanityCheck;

  $m->displayCounts() if ($dispCounts);
}



sub usage {
  print "Usage: $0 [options] [filename]\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default $_->{DEFAULT})" : "";
  }
  exit(911);
}
