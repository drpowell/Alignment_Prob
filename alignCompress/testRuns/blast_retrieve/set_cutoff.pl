#!/usr/bin/perl -w

use strict;
use Data::Dumper;

my(@Options,$DEBUG, $cutoff, $numhits);
setOptions();

(defined($cutoff) && defined($numhits)) && die "Set ONE of 'cutoff' and 'numhits";

(!defined($cutoff) && !defined($numhits)) && die "Set ONE of 'cutoff' and 'numhits";

my $hits = eval_perl_struct(join '',<>);

my $num = 0;

my $res;
for my $s (@$hits) {
  my $a = { %$s };      # copy it.
  delete $a->{SEQS};
  for my $seqs (@{$s->{SEQS}}) {
    $num++;
    ($seqs->{EVALUE} =~ /^e/) && ($seqs->{EVALUE} = "1" . $seqs->{EVALUE});
    next if (defined($cutoff) && $seqs->{EVALUE} > $cutoff);
    push(@$res, { %$a, %$seqs });
  }
}
$hits = $res;           # $hits is now flattened.

$hits = [sort { $a->{EVALUE} <=> $b->{EVALUE} } @$hits];

if (defined($numhits) && $numhits < @$hits) {
  $hits = [ @$hits[0 .. $numhits-1] ];
}

printf STDERR "Originally %d seqs.  Now %d seqs.  E-val range %g to %g\n",
              $num, scalar @$hits, $hits->[0]{EVALUE},  $hits->[-1]{EVALUE};

#map { printf "%g\n", $_->{EVALUE} } @$hits;
print Dumper($hits);


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
	      {OPT=>"cutoff=s",    VAR=>\$cutoff, DESC=>"Evalue cutoff"},
              {OPT=>"numhits=i",   VAR=>\$numhits, DESC=>"Number of hit sequences to return"},
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



