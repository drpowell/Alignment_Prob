#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use IPC::Open3;
use IO::Handle;
use IO::Select;
use IO::File;

my $timeProg = '/usr/bin/time';
my $compProg = 'java -Xmx512m -cp ../../..:../../../hb15.zip alignCompress/AlignCompress';
my $compProgOpt = ' --markov=1 --verbose=1 --maxIterations=20 --linear=true --local=true --sum=true';
my $prssProg = '../prss33 -b 200 -n -q';



my(@Options,$DEBUG,$seqFile,$qseqFile);
setOptions();

(!defined($seqFile) || !defined($qseqFile)) && die "You must specify a query seq file AND a hits seq file.";

$|=1;        # Flush stdout

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

my $str  = "";
$str .= `hostname`."\n";
$str .= `uname -a`."\n";
$str .= `free`."\n";
$str .= `date`."\n";
$str .= "pid=$$\n";
$str .= "\nProgs to use:\n$compProg$compProgOpt\n$prssProg\n";
$str .= sprintf "\nquery seq file=$qseqFile.  len=%d\n",length($qSeq);
$str .= "\nhit seq file=$seqFile.\n";

$str =~ s/^/#/gm;
print $str;

my @to_delete;

my $log = new IO::File "> outLog.$$";
#my $log = new IO::File "> /dev/null";
(defined $log) || die "Can't open output log";
$log->autoflush(1);

my $numUpto = 1;
for my $hitSeq (@$hits) {
  my $str1 = $qSeq;
  my $str2 = $hitSeq->{H_SEQ};

  printf $log "Doing %d of %d\n", $numUpto, scalar @$hits;
  print $log "s1=$str1\ns2=$str2\n";

  my($r, $rTime1, $uTime1, $sTime1, $swaps1) = 
    runProg($compProg . $compProgOpt, $str1, $str2);

  $hitSeq->{ALIGNCOMPRESS} = {r => $r->[-1], UTIME => $uTime1};

  my($prob,$score,$expect,$num, $rTime2, $uTime2, $sTime2, $swaps2) = 
      runFastaProg($prssProg, $str1, $str2);

  $hitSeq->{PRSS} = {p => $prob, SCORE=>$score, EXPECT=>$expect, RUNS=>$num, UTIME => $uTime2};

  printf "%s\n", Data::Dumper->new([$hitSeq])->Indent(0)->Dump();
  printf $log "\n%s\n\n",Dumper($hitSeq);
  $numUpto++;
}

sub runProg {
  my($prog, $str1, $str2) = @_;

  open3(\*WTR,\*RDR,\*ERR, "$timeProg $prog") ||
    die "Can't run $prog";
  print WTR "$str1\n";
  print WTR "$str2\n";

  my $rdr = new IO::Handle;
  my $err = new IO::Handle;
  $rdr->fdopen(fileno(RDR), "r");
  $err->fdopen(fileno(ERR), "r");

  my $s = new IO::Select;
  $s->add($rdr, $err);

  my(@odds_ratio);
  my($rTime,$uTime,$sTime) = (-1,-1,-1);
  my($swaps)               = (-1);

  while(my @ready = $s->can_read) {
    foreach my $fh (@ready) {
      ($fh == $rdr) && do {
        (!defined($_ = <$rdr>)) && do {$s->remove($rdr); next;};

#        printf $log "STDOUT: %s",$_;
        if (/log odds ratio = (\S+)/) { push(@odds_ratio,$1);}
        next;
      };
      ($fh == $err) && do {
        (!defined($_ = <$err>)) && do {$s->remove($err); next};

	print "STDERR: $_" if (/^NON-CONVERGENCE/);
        printf $log "STDERR:%s: %s",time(),$_;
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

sub runFastaProg {
  my($prog, $str1, $str2) = @_;

  my($f1,$f2) = ("in.$$.seq1","in.$$.seq2");

  push(@to_delete, $f1,$f2);
  $SIG{QUIT} = \&done;
  $SIG{INT} = \&done;

  open(F,"> $f1") or die "Can't create tmp file";
  print F ">GEN_IN1\n$str1\n";
  close(F);
  open(F,"> $f2") or die "Can't create tmp file";
  print F ">GEN_IN2\n$str2\n";
  close(F);

  open3(\*WTR,\*RDR,\*ERR, "$timeProg $prog $f1 $f2") ||
    die "Can't run $prog";
#  print WTR "$str1\n";
#  print WTR "$str2\n";

  my $rdr = new IO::Handle;
  my $err = new IO::Handle;
  $rdr->fdopen(fileno(RDR), "r");
  $err->fdopen(fileno(ERR), "r");

  my $s = new IO::Select;
  $s->add($rdr, $err);

  my($score,$prob,$num,$expect);
  my($rTime,$uTime,$sTime) = (-1,-1,-1);
  my($swaps)               = (-1);

  while(my @ready = $s->can_read) {
    foreach my $fh (@ready) {
      ($fh == $rdr) && do {
        (!defined($_ = <$rdr>)) && do {$s->remove($rdr); next;};

#        printf $log "STDOUT: %s",$_;
	if (/unshuffled s-w score: (\d+).*p\(\1\) < (\S+)/) {
	  $score = $1;
	  $prob  = $2;
	}
	if (/For (\d+) sequences, a score.*is expected\s+(\S+)\s+times/) {
	  $num = $1;
	  $expect = $2;
	}
        next;
      };
      ($fh == $err) && do {
        (!defined($_ = <$err>)) && do {$s->remove($err); next};

        printf $log "STDERR:%s: %s",time(),$_;
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

  (!defined($prob))  && (die "Unable to find 'prob'!\n");
  (!defined($score)) && (die "Didn't get score\n");
  (!defined($expect))&& (die "Didn't get expect\n");
  (!defined($num))   && (die "Didn't get num\n");

  unlink(@to_delete);
  pop(@to_delete);
  pop(@to_delete);
  $SIG{QUIT} = 'DEFAULT';
  $SIG{INT}  = 'DEFAULT';


  return ($prob, $score, $expect, $num, $rTime, $uTime, $sTime, $swaps);
}

sub done {
  unlink(@to_delete);
  exit;
}

sub rand_length {
  my($mid, $range) = @_;
  return $mid + int(rand($range*2+1))-$range;
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



