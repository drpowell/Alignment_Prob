#!/usr/bin/perl -w

use strict;

#use Markov_gen;
use Multi_model;

my $timeProg = '/usr/bin/time';
my $compProg = 'java -Xmx512m -cp ../..:../../hb15.zip alignCompress/AlignCompress';
my $compProgOpt = ' --markov=1 --verbose=2 --maxIterations=20 --linear=true --local=true --sum=true';
my $prssProg = './prss33 -b 200 -n -q';

use Data::Dumper;
use IPC::Open3;
use IO::Handle;
use IO::Select;
use IO::File;

$|=1;        # Flush stdout

my $log = new IO::File "> popLog2.$$";
#my $log = new IO::File "> /dev/null";
(defined $log) || die "Can't open output log";
$log->autoflush(1);

my $model = new Multi_model();

my $numArchetypes = 10;
my $numEachMutations = 5;
my @numMutation = (30, 40, 50, 60, 80);
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
$str .= "\nProgs to use:\n$compProg$compProgOpt\n$prssProg\n\n";

$str .= "Produce $numArchetypes sequences of the form gen($l1_s+-$l_range) . sub_seq($l_sub+-$l_sub_range)) . gen($l1_e+-$l_range)\n";
$str .= "From these children will be produced of the form gen($l2_s+-$l_range). mutate_sub_seq(numMutate).gen($l2_e+-$l_range)\n";
$str .= "Repeat each mutation rate $numEachMutations. Num mutations = (@numMutation)\n\n";
#$str .= "Note the model has biased 1st order stats, _but_ uniform 0 order stats\n";
$str .= "model: " . $model->as_string();
#$str .= "model2: " . $model2->as_string();

$str =~ s/^/#/gm;
print $str;

my @archetypes;
my @population;

for my $i (0 .. $numArchetypes-1) {
  my $m = $model;
  my $subseq = $m->gen_sequence(rand_length($l_sub, $l_sub_range));
  my $s1 = $m->gen_sequence(rand_length($l1_s, $l_range));
  my $e1 = $m->gen_sequence(rand_length($l1_e, $l_range));
  my $str1 = $s1 . $subseq . $e1;
  $archetypes[$i] = {SEQ => $str1, S_LEN=>length($s1), E_LEN=>length($e1), SUB_LEN=>length($subseq)};
  for my $numMutations (@numMutation) {
    for my $j (1 .. $numEachMutations) {
      my $s2 = $m->gen_sequence(rand_length($l2_s, $l_range));
      my $e2 = $m->gen_sequence(rand_length($l2_e, $l_range));
      my $subseq2 = $m->mutate($subseq, $numMutations);
      my $str2 = $s2 . $subseq2 . $e2;
      push( @population, {SEQ => $str2, PARENT=>$i, MUTATES=>$numMutations,
			  S_LEN=>length($s2), E_LEN=>length($e2), SUB_LEN=>length($subseq2)});
    }
  }
}

{
  local $Data::Dumper::Indent=3;
  print $log (Data::Dumper->Dump([\@archetypes],['$archetype']));
  print $log (Data::Dumper->Dump([\@population],['$population']));
  print $log "\n\n";
}

my @to_delete;

my $numToRun = num_processors();
my $numRunning = 0;

$SIG{CHLD} = sub { my $pid=wait;
                   if ($pid>0) {
  		     #print STDERR "Child $pid finished\n";
		     $numRunning--;
                   };
                 };

for my $i (0 .. $numArchetypes-1) {
  my $str1 = $archetypes[$i]{SEQ};
  for my $j (0 .. $#population) {
    my $str2 = $population[$j]{SEQ};

    while ($numRunning >= $numToRun) { # Wait until someone exits.
      sleep(1);
    }

    $numRunning++;
    my $pid = fork();
    die "Fork failed!" if (!defined($pid));

    next if ($pid);		# I'm the parent continue;

    # I'm the child;
    #print STDERR "Running process $$\n";
    $SIG{CHLD} = 'IGNORE';

    for my $blend (qw(false true)) {
      my($r, $a_len, $m_len, $d_len, $params, $rTime1, $uTime1, $sTime1, $swaps1) = 
	runProg($compProg . $compProgOpt . " --blend=$blend", $str1, $str2);

      printf("AlignCompress (blend=$blend): s1=%d s2=%d parent=%d mutates=%d r=%f (%f) al=%f (%f) ml=%f (%f) dl=%f (%f) uTime=%f params:%s\n",
	     $i, $j,
	     $population[$j]{PARENT},
	     $population[$j]{MUTATES},
	     $r->[-1], $r->[0],
	     $a_len->[-1], $a_len->[0],
	     $m_len->[-1], $m_len->[0],
	     $d_len->[-1], $d_len->[0],
	     $uTime1,
	     join " ", map { "$_=$params->{$_} " } sort keys %$params,
	    );
    }

    my($prob,$score,$expect,$num, $rTime2, $uTime2, $sTime2, $swaps2) = 
      runFastaProg($prssProg, $str1, $str2);
    printf("PRSS: s1=%d s2=%d parent=%d mutates=%d p=$prob s=$score expect=$expect " .
	   "runs=$num uTime=$uTime2\n", 
	   $i, $j, 
	   $population[$j]{PARENT},
	   $population[$j]{MUTATES});

    printf $log "\nDONE\n\n";

    exit;			# The child is done
  }
}

while ($numRunning>0) {
  sleep(1);
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

  my(@odds_ratio, @total_len, @model_len, @data_len);
  my(%params);
  my($rTime,$uTime,$sTime) = (-1,-1,-1);
  my($swaps)               = (-1);

  while(my @ready = $s->can_read) {
    foreach my $fh (@ready) {
      ($fh == $rdr) && do {
        (!defined($_ = <$rdr>)) && do {$s->remove($rdr); next;};

#        printf $log "STDOUT: %s",$_;
        if (/log odds ratio = (\S+)/)  { push(@odds_ratio,$1);}
	if (/Mutual Encoding = ([.\d]+).*model=([.\d]+) data=([.\d]+)/) {
	  push(@total_len, $1);
	  push(@model_len, $2);
	  push(@data_len,  $3);
	}

	if (/^(diag_fromD   |
	       start_fromD  |
	       diag_fromI   |
	       start_fromI  |
	       cont_fromI   |
	       match_cost   |
	       change_cost)  =  (\S+)/x) {
	  $params{$1} = $2;
	}

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

  return (\@odds_ratio, \@total_len, \@model_len, \@data_len, \%params, $rTime, $uTime, $sTime, $swaps);
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

sub num_processors {
  open(F, "< /proc/cpuinfo") or do { print STDERR "Can't read /proc/cpuinfo\n"; return 1; };
  my(@l) = <F>;
  my $num = grep (/^processor/, @l);
  if ($num<=0) {
    print STDERR "Bad number of processors: $num\n";
    return 1;
  }
  return $num;
}

