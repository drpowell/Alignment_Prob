#!/usr/bin/perl -w

use strict;

my($myLoc) = ($0 =~ m{^(.*)/});
$myLoc ||= '.';

my $fname = shift;
my $MMorder="unknown";
my $entropy="unknown";
my $pMM="unknown";
my $uni0 = 0;
my $haveSW = 0;

open (F, "< $fname") or die "Unable to read $fname";
while (<F>) {
  if (/^\#Model order=(\S+) entropy=(\S+)/) {
    $MMorder=$1;
    $entropy=$2;
  }

  if (/^\#java.*--markov=(\S+)/) {
    $pMM=$1;
  }

  if (/^\#Note the model has biased 1st order stats, _but_ uniform 0 order stats$/) {
    $uni0 = 1;
  }

  if (/^AlignCompress \(SW\)/) {
#    $haveSW = 1;
  }
}

my $cmd = "";

my $toFile = 0;

if (defined($ARGV[0])) {
  $toFile = 1;
  my $file = shift;
  die "$file doesn't have .eps extension" if (!($file =~ /\.eps$/));
  $cmd .= "set terminal postscript eps\n";
  $cmd .= "set output '$file'\n";
}

$cmd .= "set key left\n";
$cmd .= "set logscale y\n";
$cmd .= sprintf("set title 'ROC for MM order=$MMorder entropy=$entropy%s'\n",
		$uni0 ? " (uniform 0 order stats)" : "");
$cmd .= "set xlabel 'Coverage'\n";
$cmd .= "set ylabel 'Errors'\n";


$cmd .= "plot '< $myLoc/extractPopData.pl prss $fname' title 'PRSS p-value' with lines";
$cmd .= ",'< $myLoc/extractPopData.pl prss2 $fname' title 'PRSS raw-score' with lines";
$cmd .= ",'< $myLoc/extractPopData.pl al_one $fname' title 'Optimal Alignment -markov=$pMM' with lines";
$cmd .= ",'< $myLoc/extractPopData.pl al_all $fname' title 'Average Alignment -markov=$pMM' with lines";
$cmd .= ",'< $myLoc/extractPopData.pl al_sw $fname' title 'Optimal SW alignment' with lines" if ($haveSW);
$cmd .= "\n";

if ($toFile) {
  $cmd .= "quit\n";
  open(F,"| gnuplot") or (die "Can't run gnuplot");
  print F $cmd;
  close(F);
} else {
  $cmd .= "pause -1\n";
  print $cmd;
}
