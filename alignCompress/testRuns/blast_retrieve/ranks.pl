#!/usr/bin/perl -w

use strict;
use Data::Dumper;

my(@Options,$DEBUG,$seqFile, $qseqFile);
setOptions();

my $hits = eval_perl_file($seqFile);

my $qSeq = readFASTA($qseqFile);
sub readFASTA {
  my($f) = @_;
  open(F, "<$f") or die "Can't read $f";
  my $str = join"",<F>;
  close(F);
  $str =~ /\A>/ or die "$f not in FASTA format";
  $str =~ s/\A>.*?\n//;     # Nuke name line
  $str =~ s/\s//g;          # Nuke all newlines and whitespace
  return $str;
}


my $blast  = [sort { $a->{EVALUE} <=> $b->{EVALUE} } @$hits];
my $Malign = [sort { $b->{ALIGNCOMPRESS}{r} <=> $a->{ALIGNCOMPRESS}{r} } @$hits];
my $prss   = [sort { $a->{PRSS}{p} <=> $b->{PRSS}{p} } @$hits];

for my $i (0 .. $#$hits) {
  $blast->[$i]{RANKS}{BLAST} = $i;
  $Malign->[$i]{RANKS}{ALIGNCOMPRESS} = $i;
  $prss->[$i]{RANKS}{PRSS} = $i;
}

print "<TABLE BORDER=1 CELLPADDING=2 CELLSPACING=2>\n";
print "<TR><TH><TH>Name<TH>Description<TH>Blast rank<TH>PRSS rank<TH>M-alignment rank\n";
for my $i (0 .. $#$hits) {
  my $s = $hits->[$i];
  printf "<TR><TD>%3d<TD><A HREF='%s'>%s</A><TD> %s  <TD><B>%d</B> (%g)<TD><B>%d</B> (%g) <TD><B>%d</B> (%g)\n",
    $i,
    sprintf("align.pl.cgi?s1=%s&s2=%s",$qSeq,$s->{H_SEQ}),
    $s->{NAME}, $s->{DESCR},
    $s->{RANKS}{BLAST}, $s->{EVALUE},
    $s->{RANKS}{PRSS},  $s->{PRSS}{p},
    $s->{RANKS}{ALIGNCOMPRESS}, $s->{ALIGNCOMPRESS}{r};
}
print "</TABLE>\n";

#print Dumper($hits);


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


#----------------------------------------------------------------------
# option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
              {OPT=>"help",      VAR=>\&usage,       DESC=>"This help"},
              {OPT=>"debug=i",   VAR=>\$DEBUG,  DEFAULT=>1,
               DESC=>"Debug level"},
              {OPT=>"seqs=s",    VAR=>\$seqFile,  DESC=>"File containing sequences (must be a perl structure)"},
              {OPT=>"qseq=s",    VAR=>\$qseqFile,  DESC=>"File containing query sequence"},
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



