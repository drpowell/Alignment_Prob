#!/usr/bin/perl -w

use strict;

my($myLoc) = ($0 =~ m{^(.*)/});
$myLoc ||= '.';

my $fname = shift;
my $MMorder="unknown";
my $entropy="unknown";
my $uni0 = 0;

open (F, "< $fname") or die "Unable to read $fname";
while (<F>) {
  if (/^\#Model order=(\S+) entropy=(\S+)/) {
    $MMorder=$1;
    $entropy=$2;
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
  $cmd .= "set terminal postscript eps\n";
  $cmd .= "set output '$file'\n";
}

$cmd .= "set key right\n";
$cmd .= sprintf("set title 'Log-odds for MM order=$MMorder entropy=$entropy%s'\n",
		$uni0 ? " (uniform 0 order stats)" : "");
$cmd .= "set xlabel 'No. mutates'\n";
$cmd .= "set ylabel 'bits'\n";


$cmd .= "plot '< $myLoc/extractData.pl prss $fname' thru -log(x)/log(2) title 'PRSS -log(p)' with linespoints";
$cmd .= ",'< $myLoc/extractData.pl al_one $fname' title 'Optimal Alignment' with linespoints";
$cmd .= ",'< $myLoc/extractData.pl al_all $fname' title 'Average Alignment' with linespoints";
$cmd .= ",0 notitle";
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
