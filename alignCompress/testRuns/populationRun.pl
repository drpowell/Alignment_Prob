#!/usr/bin/perl -w

use strict;

use Markov_gen;

my $timeProg = '/usr/bin/time';
my $compProg = 'java -Xmx512m -cp ../..:../../hb15.zip alignCompress/AlignCompress';
my $compProgOpt = ' --markov=-1 --verbose=1 --maxIterations=20 --linear=true --local=true';
my $prssProg = './prss33 -b 200 -n -q';

use Data::Dumper;
use IPC::Open3;
use IO::Handle;
use IO::Select;
use IO::File;

$|=1;        # Flush stdout

my $log = new IO::File "> popLog.$$";
#my $log = new IO::File "> /dev/null";
(defined $log) || die "Can't open output log";
$log->autoflush(1);

my $model = new Markov_gen(-1, [qw(a t g c)]);
#$model->model_power(2);

my $numArchetypes = 10;
my $numEachMutations = 4;
my @numMutation = (10, 30, 50, 80, 100);

my $str  = "";
$str .= `hostname`."\n";
$str .= `uname -a`."\n";
$str .= `free`."\n";
$str .= `date`."\n";
$str .= "pid=$$\n";
$str .= "\nProgs to use:\n$compProg$compProgOpt\n$prssProg\n\n";

$str .= "Produce $numArchetypes sequences of the form gen(200) . sub_seq(100) . gen(100)\n";
$str .= "From these children will be produced of the form gen(100). mutate_sub_seq(numMutate).gen(200)\n";
$str .= "Repeat each mutation rate $numEachMutations. Num mutations = (@numMutation)\n\n";
$str .= $model->as_string();

$str =~ s/^/#/gm;
print $str;

my @archetypes;
my @population;

for my $i (0 .. $numArchetypes-1) {
  my $subseq = $model->gen_sequence(100);
  my $str1 = $model->gen_sequence(200) . $subseq . $model->gen_sequence(100);
  $archetypes[$i] = $str1;
  for my $numMutations (@numMutation) {
    for my $j (1 .. $numEachMutations) {
      my $str2 = $model->gen_sequence(100) . $model->mutate($subseq, $numMutations) .
	$model->gen_sequence(200);
      push( @population, {SEQ => $str2, PARENT=>$i, MUTATES=>$numMutations});
    }
  }
}

{
  local $Data::Dumper::Indent=3;
  print $log (Data::Dumper->Dump([\@archetypes],['@archetype']));
  print $log (Data::Dumper->Dump([\@population],['@population']));
  print $log "\n\n";
}

my @to_delete;

for my $i (0 .. $numArchetypes-1) {
  my $str1 = $archetypes[$i];
  for my $j (0 .. $#population) {
    my $str2 = $population[$j]{SEQ};

    print $log "s1=$str1\ns2=$str2\n";

    for my $sum (qw(true false)) {
      my($r, $rTime1, $uTime1, $sTime1, $swaps1) = 
	runProg($compProg . $compProgOpt . " --sum=$sum", $str1, $str2);

      printf("AlignCompress (sum=$sum): s1=%d s2=%d parent=%d mutates=%d r=%f uTime=%f\n",
	     $i, $j,
	     $population[$j]{PARENT},
	     $population[$j]{MUTATES},
	     $r->[-1],$uTime1);
    }

    my($r, $rTime1, $uTime1, $sTime1, $swaps1) = 
	runProg($compProg . $compProgOpt . " --sum=true --doSW=true", $str1, $str2);

      printf("AlignCompress (SW): s1=%d s2=%d parent=%d mutates=%d r=%f uTime=%f\n",
	     $i, $j,
	     $population[$j]{PARENT},
	     $population[$j]{MUTATES},
	     $r->[-1],$uTime1);

    my($prob,$score,$expect,$num, $rTime2, $uTime2, $sTime2, $swaps2) = 
      runFastaProg($prssProg, $str1, $str2);
    printf("PRSS: s1=%d s2=%d parent=%d mutates=%d p=$prob s=$score expect=$expect " .
	   "runs=$num uTime=$uTime2\n", 
	   $i, $j, 
	   $population[$j]{PARENT},
	   $population[$j]{MUTATES});

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
        printf $log "STDERR: %s",$_;
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

        printf $log "STDERR: %s",$_;
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