#!/usr/bin/perl -w

use strict;

#use Markov_gen;
use Multi_model;

my $timeProg = '/usr/bin/time';
my $compProg = 'java -Xmx512m -cp ../..:../../hb15.zip alignCompress/AlignCompress';
my $compProgOpt = ' --markov=0 --verbose=1 --maxIterations=20 --linear=true --local=true';
my $prssProg = './prss33 -b 200 -n -q';

use IPC::Open3;
use IO::Handle;
use IO::Select;
use IO::File;

$|=1;        # Flush stdout

my $log = new IO::File "> outLog.$$";
#my $log = new IO::File "> /dev/null";
(defined $log) || die "Can't open output log";
$log->autoflush(1);

my $model = new Multi_model();
#my $model = new Markov_gen(-1, [qw(a t g c)]);
#$model->model_power(2);
#$model->makeUniModel(2);


my ($l_sub, $l_sub_range) = (120, 30);
my ($l1_s, $l1_e) = (50, 100);
my ($l2_s, $l2_e) = (100, 50);
my ($l_range) = 50;


my $str  = "";
$str .= `hostname`."\n";
$str .= `uname -a`."\n";
$str .= `free`."\n";
$str .= `date`."\n";
$str .= "pid=$$\n";
$str .= "\nProgs to use:\n$compProg$compProgOpt\n$prssProg\n";
#$str .= "\nSequences are unrelated both of length exactly 400 characters\n\n";
$str .= "\nSequences are related.\n";
$str .= "subseq1 = gen($l_sub +- $l_sub_range)\n";
$str .= "subseq2 = mutate(subseq1, n_times)\n";
$str .= "s1 = gen($l1_s+-$l_range) . subseq1 . gen($l1_e+-$l_range)\n";
$str .= "s2 = gen($l2_s+-$l_range) . subseq2 . gen($l2_e+-$l_range)\n\n";
#$str .= "Note the model has biased 1st order stats, _but_ uniform 0 order stats\n";
$str .= $model->as_string();

$str =~ s/^/#/gm;
print $str;

my $numRuns = 10;

my @to_delete;

for my $numMutations (0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 300) {
  for my $runNum (1..$numRuns) {
    my $subseq = $model->gen_sequence(rand_length($l_sub, $l_sub_range));
    my $str1 = $model->gen_sequence(rand_length($l1_s, $l_range)) . $subseq . $model->gen_sequence(rand_length($l1_e, $l_range));
    my $str2 = $model->gen_sequence(rand_length($l2_s, $l_range)) . $model->mutate($subseq, $numMutations) . $model->gen_sequence(rand_length($l2_e, $l_range));
#    my $str1 = $model->gen_sequence(400);
#    my $str2 = $model->gen_sequence(400);
#    my $numMutations = -1; # Use -1 to denote unrelated sequences

    print $log "s1=$str1\ns2=$str2\n";

    for my $blend (qw(true false)) {
      for my $sum (qw(true false)) {
	my($r, $rTime1, $uTime1, $sTime1, $swaps1) = 
	  runProg($compProg . $compProgOpt . " --sum=$sum --blend=$blend", $str1, $str2);

	printf "AlignCompress (blend=$blend sum=$sum): mutates=$numMutations r=%f uTime=%f\n",$r->[-1],$uTime1;
      }
    }

    my($prob,$score,$expect,$num, $rTime2, $uTime2, $sTime2, $swaps2) = 
      runFastaProg($prssProg, $str1, $str2);
    print "PRSS:  mutates=$numMutations p=$prob s=$score expect=$expect runs=$num uTime=$uTime2\n";

    printf $log "\nDONE\n\n";
  }
}


sub runProg {
  my($prog, $str1, $str2) = @_;

  open3(undef,\*RDR,\*ERR, "$timeProg $prog $str1 $str2") ||
    die "Can't run $prog";
#  print WTR "$str1\n";
#  print WTR "$str2\n";

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

  open3(undef,\*RDR,\*ERR, "$timeProg $prog $f1 $f2") ||
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
