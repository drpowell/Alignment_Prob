#!/usr/bin/perl -w

use strict;
use Data::Dumper;

my(@Options,$DEBUG, $file1, $file2, $order);
setOptions();

(!defined($file1) || !defined($file2)) && do {
  print "You must specify file1 AND file2\n";
  usage();
};

open(F1, "<$file1") or die "Can't open '$file1'";
open(F2, "<$file2") or die "Can't open '$file2'";

while (<F1>) {
  if (/^#/ || /^\s*$/) {
    print;
    next;
  }
  my $s1 = eval_perl_str($_);
  my $s2;
  while (1) {
    my $line = <F2>;
    last if !defined($line);
    unless ($line=~/^#/ || $line=~/^\s*$/) {
      $s2 = eval_perl_str($line);
      last;
    }
  }

  $s2 = $s1 if (!defined($s2));

  # Okay merge $s1 and $s2
  # First check they are for the same seq.
  die "Name mismatch!" if ($s1->{NAME} ne $s2->{NAME});
  die "Start pos mismatch!" if ($s1->{H_START} != $s2->{H_START});
  die "End pos mismatch!" if ($s1->{H_END} != $s2->{H_END});

  # For now, just put {ALIGNCOMPRESS2} into $s1
  die "'$file1' already contains 'other' results" if exists($s1->{'ALIGNCOMPRESS' . $order});
  $s1->{'ALIGNCOMPRESS' . $order} = $s2->{ALIGNCOMPRESS};
  printf "%s\n", Data::Dumper->new([$s1])->Indent(0)->Dump();
}


sub eval_perl_str {
  my($str) = @_;
  my($VAR1);
  my $res = eval($str);
  if ($@) { die "Error evaling : $@"; };
  return $res;
}


#----------------------------------------------------------------------
# option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
              {OPT=>"help",      VAR=>\&usage,       DESC=>"This help"},
              {OPT=>"debug=i",   VAR=>\$DEBUG,  DEFAULT=>1,
               DESC=>"Debug level"},
              {OPT=>"f1=s",    VAR=>\$file1, DESC=>"First file"},
              {OPT=>"f2=s",   VAR=>\$file2, DESC=>"Second file"},
	      {OPT=>'order=i',  VAR=>\$order, DEFAULT=>2, 
	       DESC=>"Integer to specify what to use in the key for the second file"},
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


