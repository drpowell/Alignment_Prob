#!/usr/bin/perl -w

use strict;

my($myLoc) = ($0 =~ m{^(.*)/});
$myLoc ||= '.';

my $fname = shift;
my(@MMorder) = ();
my(@entropy) = ();
my(@Pchange) = ();
my $pMM="unknown";
my $uni0 = 0;
my $haveSW = 0;

open (F, "< $fname") or die "Unable to read $fname";
while (<F>) {
  if (/Model order=(\S+) entropy=(\S+) Pchange=(\S+)/) {
    push(@MMorder, $1);
    push(@entropy, $2);
    push(@Pchange, $3);
  }

  if (/^\#java.*--markov=(\S+)/) {
    $pMM=$1;
  }

  if (/^\#Note the model has biased 1st order stats, _but_ uniform 0 order stats$/) {
    $uni0 = 1;
  }
}

my $cmd = "";

my $toFile = 0;

if (defined($ARGV[0])) {
  $toFile = 1;
  my $file = shift;
  die "$file doesn't have .eps extension" if (!($file =~ /\.eps$/));
#  $cmd .= "set terminal postscript eps color\n";
  $cmd .= "set terminal postscript eps\n";
  $cmd .= "set output '$file'\n";
}

$cmd .= "set key left\n";
$cmd .= "set logscale y\n";
$cmd .= "set data style lines\n";
if (@MMorder == 1) {
  $cmd .= "set title \"ROC for MM order=$MMorder[0] entropy=$entropy[0] Pchange=$Pchange[0]";
} else {
  $cmd .= "set title \"ROC for blend of MM order=$MMorder[0] entropy=$entropy[0] Pchange=$Pchange[0]\\n".
                          "and MM order=$MMorder[1] entropy=$entropy[1] Pchange=$Pchange[1]";
}
$cmd .= ($uni0 ? " (uniform 0 order stats)" : "") . "\"\n";
$cmd .= "set xlabel 'Coverage'\n";
$cmd .= "set ylabel 'Errors'\n";


$cmd .= "plot '< $myLoc/extractPopData.pl blend prss $fname' title 'PRSS p-value'";
$cmd .= ",'< $myLoc/extractPopData.pl blend prss2 $fname' title 'PRSS raw-score'";
#$cmd .= ",'< $myLoc/extractPopData.pl al_one $fname' title 'Optimal Alignment -markov=$pMM'";
$cmd .= ",'< $myLoc/extractPopData.pl blend al_one $fname' title 'Average Alignment -markov=$pMM'";
$cmd .= ",'< $myLoc/extractPopData.pl blend al_all $fname' title 'Average Alignment -blendModel'";
#$cmd .= ",'< $myLoc/extractPopData.pl al_sw $fname' title 'Optimal SW alignment'" if ($haveSW);
$cmd .= "\n";

if ($toFile) {
  $cmd .= "quit\n";
  open(F,"| gnuplot") or (die "Can't run gnuplot");
  print F $cmd;
  close(F);
} else {
  open(F, "| gnuplot -persist -") or (die "Can't run gnuplot");
  { my $oldfh = select F; $|=1; select $oldfh; }
  print F $cmd;
  close(F);
#  print $cmd;
}
