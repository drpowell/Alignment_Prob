#!/usr/bin/perl -w

use strict;

my($myLoc) = ($0 =~ m{^(.*)/});
$myLoc ||= '.';

my $fname = shift;
my @MMorder="unknown";
my @entropy="unknown";
my @Pchange;
my $uni0 = 0;
my $blend = 0;

open (F, "< $fname") or die "Unable to read $fname";
while (<F>) {
  if (/^\#Model order=(\S+) entropy=(\S+)/) {
    $MMorder[0]=$1;
    $entropy[0]=$2;
  }

  if (/^\#Blend Model:/) {
    $blend=1;
  }

  if ($blend && /Model(\d): Model order=(\S+) entropy=(\S+) Pchange=(\S+)/) {
    $MMorder[$1-1] = $2;
    $entropy[$1-1] = $3;
    $Pchange[$1-1] = $4;
  }

  if (/^\#Note the model has biased 1st order stats, _but_ uniform 0 order stats$/) {
    $uni0 = 1;
  }

}

my $cmd = "";

map { $_ = 'uniform' if ($_==-1) } @MMorder;

my $toFile = 0;

if (defined($ARGV[0])) {
  $toFile = 1;
  my $file = shift;
  die "$file doesn't have .eps extension" if (!($file =~ /\.eps$/));
  $cmd .= "set terminal postscript eps\n";
  $cmd .= "set output '$file'\n";
}

$cmd .= "set key right\n";
$cmd .= sprintf('set title "Log-odds for MM order=%d entropy=%.4f %s%s%s"'."\n",
		$MMorder[0], $entropy[0],
		(@Pchange ? "Pchange=$Pchange[0]" : ''),
		($blend ? "\\nand MM order=$MMorder[1] entropy=$entropy[1] Pchange=$Pchange[1]" : ""),
		$uni0 ? " (uniform 0 order stats)" : "");
$cmd .= "set xlabel 'No. mutates'\n";
$cmd .= "set ylabel 'bits'\n";

$cmd .= "pf(x) = -log(x)/log(2)\n";

$cmd .= "plot '< $myLoc/extractData.pl prss $fname' using 1:(pf(\$2)) title 'PRSS -log(p)' with linespoints";
#$cmd .= ",'< $myLoc/extractData.pl prss $fname' using 1:(pf(\$2)):(pf(\$3)):(pf(\$4)) notitle with errorbars ";
#$cmd .= ",'< $myLoc/extractData.pl al_one $fname' title 'Optimal Alignment' with linespoints";
#$cmd .= ",'< $myLoc/extractData.pl al_one $fname' notitle with errorbars";
$cmd .= ",'< $myLoc/extractData.pl al_all $fname' title 'Average Alignment' with linespoints";
#$cmd .= ",'< $myLoc/extractData.pl al_all $fname' notitle with errorbars";
#$cmd .= ",'< $myLoc/extractData.pl al_one_blend $fname' title 'Optimal Alignment - blend model' with linespoints" if ($blend);
#$cmd .= ",'< $myLoc/extractData.pl al_one_blend $fname' notitle with errorbars" if ($blend);
$cmd .= ",'< $myLoc/extractData.pl al_all_blend $fname' title 'Average Alignment - blend model' with linespoints" if ($blend);
#$cmd .= ",'< $myLoc/extractData.pl al_all_blend $fname' notitle with errorbars" if ($blend);
$cmd .= ",0 notitle";
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

