#!/usr/bin/perl -w

use strict;
use IPC::Open3;
use IO::Handle;
use IO::Select;
use IO::File;


my $str = join(" ", @ARGV);

my $timeProg = '/usr/bin/time';
my $compProg = 'java -Xmx512m -cp ../..:../../hb15.zip alignCompress/AlignCompress';

for my $linear qw(false true) {
  for my $sum qw(false true) {
    for my $local qw(false true) {
      my($r1) = runProg($compProg . " --markov=0 --iterations=9" 
                                  . " --linear=$linear"
                                  . " --sum=$sum"
                                  . " --local=$local", $str);

      printf("local=%-5s sum=%-5s linear=%-5s %s\n",
              $local,$sum,$linear,
	      join(" ", map { sprintf "%.4f",$_ } @$r1));
    }
  }
}


sub runProg {
  my($prog, $str) = @_;

  open3(undef,\*RDR,\*ERR, "$timeProg $prog $str") ||
    die "Can't run $prog";
#  print WTR "$str1\n";
#  print WTR "$str2\n";

  my $rdr = new IO::Handle;
  my $err = new IO::Handle;
  $rdr->fdopen(fileno(RDR), "r");
  $err->fdopen(fileno(ERR), "r");

  my $s = new IO::Select;
  $s->add($rdr, $err);

  my(@odds_ratio)          = ();
  my($rTime,$uTime,$sTime) = (-1,-1,-1);
  my($swaps)               = (-1);

  while(my @ready = $s->can_read) {
    foreach my $fh (@ready) {
      ($fh == $rdr) && do {
        (!defined($_ = <$rdr>)) && do {$s->remove($rdr); next;};

        #printf "STDOUT: %s",$_;
        if (/log odds ratio = (\S+)/) { push(@odds_ratio,$1);}
        next;
      };
      ($fh == $err) && do {
        (!defined($_ = <$err>)) && do {$s->remove($err); next};

        #printf "STDERR: %s",$_;
        if (/^([^\s]*?)user/) { $uTime = $1 };
        if (/\s([^\s]*?)system/) { $sTime = $1 };
        if (/\s([^\s]*?)elapsed/)  { $rTime = $1 };
        if (/\s([^\s]*?)swaps/) {$swaps=$1};
        next;
      };
      die "Ack bad handle to read";
    }
  }

  wait;

  (!@odds_ratio) && (die "Unable to find 'log odds ratio'!\n");

  return (\@odds_ratio, $rTime, $uTime, $sTime, $swaps);
}


