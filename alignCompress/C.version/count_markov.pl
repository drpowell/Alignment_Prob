#!/usr/bin/perl -w

use strict;
my(@Options,$DEBUG,$order,$seqFile,$qseqFile);
setOptions();

(!defined($seqFile)) && die "You must specify a hits seq file.";


my $hits = eval_perl_file($seqFile);

sub eval_perl_file {
  my($f) = @_;
  open(F, "<$f") or die "Can't read $f";
  my $str = join"",<F>;
  close(F);
  my($VAR1);
  my $res = eval($str);
  if ($@) { die "Error evaling : $@"; };
  return $res;
}

my (%c,%t);
for my $h (@$hits) {
  my $s = $h->{H_SEQ};

  for my $i (0 .. length($s)-$order-1) {
    $c{substr($s, $i, $order+1)}++;
#    $t{substr($s, $i, $order)}++;
  }
}

for my $k (sort keys %c) {
  next if $k !~ /^[atgc]*$/i;
  print "$k: $c{$k}\n";
}







#----------------------------------------------------------------------
# option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
              {OPT=>"help",      VAR=>\&usage,       DESC=>"This help"},
              {OPT=>"debug=i",   VAR=>\$DEBUG,  DEFAULT=>1,
               DESC=>"Debug level"},
              {OPT=>"seqs=s",    VAR=>\$seqFile,  DESC=>"File containing sequences (must be a perl structure)"},
              {OPT=>"order=i",   VAR=>\$order, DEFAULT=>1,
               DESC=>"Order of Markov model to fit"},

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

