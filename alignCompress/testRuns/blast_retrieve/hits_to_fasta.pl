#!/usr/bin/perl -w

use strict;
use Data::Dumper;

my(@Options,$DEBUG);
setOptions();

my $hits = eval_perl_struct(join '',<>);

my $num = 0;

my $res;
for my $s (@$hits) {
  printf ">%s %s\n", $s->{NAME}, $s->{DESCR};
  printf "%s\n\n", $s->{H_SEQ};
}


sub eval_perl_struct {
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



