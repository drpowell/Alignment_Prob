#!/usr/bin/perl -w
#
# Run the ssearch33 program that comes with fasta
#
# Written by David Powell 27/8/2002

my $ssearchProg = './ssearch33 -q';

use strict;
use IPC::Open3;
use IO::Handle;
use IO::Select;
use IO::File;

if (@ARGV != 2) {
  die "Usage: $0 file1 file2";
}


my($f1,$f2) = (shift, shift);

my($alignments) = runFastaProg($f1,$f2, $ssearchProg);

for my $i (0..$#$alignments) {
  print "ALIGNMENT:\n" . $alignments->[$i];
}

sub runFastaProg {
  my($f1, $f2, $prog) = @_;

  open3(undef,\*RDR,\*ERR, "$prog $f1 $f2") ||
    die "Can't run $prog";

  my $rdr = new IO::Handle;
  my $err = new IO::Handle;
  $rdr->fdopen(fileno(RDR), "r");
  $err->fdopen(fileno(ERR), "r");

  my $s = new IO::Select;
  $s->add($rdr, $err);

  my $buf;
  my @align;
  my $inAlign = 0;
  my $blank = 0;

  while(my @ready = $s->can_read) {
    foreach my $fh (@ready) {
      ($fh == $rdr) && do {
        (!defined($_ = <$rdr>)) && do {$s->remove($rdr); next;};

        if (/^\s*$/) {
          $blank++;

          if ($blank>=2) {
            push(@align, $buf) if ($buf && $inAlign);
            $inAlign = 0;
          }
        } else {
          $blank = 0;
        }

        if (/^>/) {
          push(@align, $buf) if $buf;
          $buf = $_;
          $inAlign = 1;
          next;
        }

        $buf .= $_ if ($inAlign);

        next;
      };
      ($fh == $err) && do {
        (!defined($_ = <$err>)) && do {$s->remove($err); next};

        printf STDERR "STDERR: %s",$_;
        next;
      };
      die "Ack bad handle to read";
    }
  }

  wait;

  return (\@align);
}

