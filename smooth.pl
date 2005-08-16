#!/usr/bin/perl -w

# D. Powell 8/8/2005
# Program to smooth the output of fuzzyLZ for plotting.
# Smoothing is done by averaging over a fixed sized window

use strict;

my(@Options,$DEBUG, $winSize, $winPos, $endFunc);
setOptions();

if (!defined($winSize)) {
  print "winSize must be set\n";
  usage();
}

my $sum;
my @window;

my @lines;

sub end_val {
  my($v) = @_;
  if ($endFunc==0) {
    return $v;
  } elsif ($endFunc==1) {
    return $sum/@window;
  } elsif ($endFunc==2) {
    return undef;
  }
}

while(<>) {
  m/^(.*?\s)?(\S+)$/ or die "Line not of expected format:\n$_";
  my($other,$v)=($1,$2);
  $other="" unless defined($other);

  push(@window, $v);
  $sum += $v;

  if (@window > $winSize) {
    $sum -= shift(@window);
  }

  if ($winPos == 0) {
    push(@lines, [$other, $v]);
    if (@lines == $winSize) {
      ($other, $v) = @{ shift(@lines) };
    }
  }
    

  my $pr_val;
  if (@window == $winSize) {
    $pr_val= $sum/@window;
  } elsif ($winPos == 1) {
    $pr_val = end_val($v);
  }

  print "$other$pr_val\n" unless (!defined($pr_val));
}

for my $line (@lines) {
  my($other, $v) = @$line;

  $sum -= shift(@window);

  my $pr_val = end_val($v);

  print "$other$pr_val\n" unless (!defined($pr_val));
}


#----------------------------------------------------------------------
# option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
            {OPT=>"help",      VAR=>\&usage,       DESC=>"This help"},
            {OPT=>"debug=i",   VAR=>\$DEBUG,  DEFAULT=>1, DESC=>"Debug level"},
            {OPT=>"winSize=i", VAR=>\$winSize,     DESC=>"Window size to average over"},
        
            {OPT=>"winPos=i",  VAR=>\$winPos, DEFAULT=>0,  DESC=>"Postion of value in smoothing window.\n    0 - Value at start of window, 1 - end"},

            {OPT=>"ends=i",    VAR=>\$endFunc,  DEFAULT=>0,  DESC=>"Functionality when smoothing window not full.\n    0 - No smoothing, 1 - average window contents, 2 - output nothing"},
             );

  #(!@ARGV) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [options]\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}

