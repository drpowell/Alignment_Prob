#!/usr/bin/perl -w

use Data::Dumper;
use strict;

use Bio::DB::GenBank;
use Bio::Tools::Run::RemoteBlast;

my $gb = new Bio::DB::GenBank();

my $end_region = 50;      # How many characters at each end to catch

  my $v = 1;
  my $prog = 'blastn';
#  my $e_val= '1e-4';
#  my $e_val= '1e-2';
  my $e_val= '10';

  my @params = ( '-prog' => $prog,
  		 '-expect' => $e_val);

  my $factory = Bio::Tools::Run::RemoteBlast->new(@params);

#  $Bio::Tools::Run::RemoteBlast::HEADER{FILTER} = '';     # Turn off low-complexity filtering

  my $results = [];
  
  $v = 1;
  my $str = Bio::SeqIO->new(-file=>'PfEMP1.bits.fa' , '-format' => 'fasta' );
  my $input = $str->next_seq();
  #  Blast a sequence against a database:
  my $r = $factory->submit_blast($input);
  print STDERR "waiting..." if( $v > 0 );
  while ( my @rids = $factory->each_rid ) {
      foreach my $rid ( @rids ) {
  	  my $rc = $factory->retrieve_blast($rid);
  	  if( !ref($rc) ) {
  	      if( $rc < 0 ) { 		
  		      $factory->remove_rid($rid);
  		  }
  	      print STDERR "." if ( $v > 0 );
  	      sleep 5;
  	  } else { 
  	      $factory->remove_rid($rid);
  	      my $result = $rc->next_result;
  	      print STDERR "db is ", $result->database_name(), "\n";
  	      my $count = 0;
  	      while( my $hit = $result->next_hit ) {		
  		  $count++;
  		  next unless ( $v > 0);
#  		  print "hit name is ", $hit->name, "\n";
#  		  print "  descr: ", $hit->description, "\n";
#  		  print "  accession:", $hit->accession, "\n";
		  push(@$results, {NAME  => $hit->name,
		                   DESCR =>  $hit->description,
				   ACCESSION => $hit->accession,
				   SEQS => [],
				   });

                  my($id) = ($hit->name =~ /\|(.*?)\|/);
		  my $seqGB;
		  while (1) {
		    $seqGB = $gb->get_Seq_by_id($id);
		    last if defined($seqGB);
		    print STDERR "FAILED to get seq id=$id     Retrying...\n";
		  }
#		  print "Retrieved seq: $seqGB\n";

  		  while( my $hsp = $hit->next_hsp ) {
#  		      printf " (%d..%d) with (%d..%d) score=%f E-val=%s\n", 
#		             $hsp->query->location->start, $hsp->query->location->end,
#		             $hsp->hit->location->start, $hsp->hit->location->end,
#		             $hsp->score, $hsp->evalue;
		      my $seq2 = lc $hsp->hit_string;
		      $seq2 =~ tr/-//d;
		      if ($hsp->hit->location->start > length($seqGB->seq())) {
		        print STDERR "BAD START POS\n";
	              }
		      my $seq1 = $seqGB->trunc($hsp->hit->location->start,$hsp->hit->location->end);
		      $seq1 = $seq1->revcom() if ($hsp->hit->location->strand < 0);
		      $seq1 = lc $seq1->seq();
		      if ($seq1 ne $seq2) {
			    print STDERR "MISMATCH on retrieved substring\n";
		      }

		      my $start = $hsp->hit->location->start - $end_region;
		      my $end   = $hsp->hit->location->end + $end_region;
		      $start = 1 if ($start < 1);
		      $end = length($seqGB->seq()) if ($end > length($seqGB->seq()));
		      my $ss = $seqGB->trunc($start, $end);
		      $ss = $ss->revcom() if ($hsp->hit->location->strand < 0);

		      push(@{$results->[-1]{SEQS}},
		                 {Q_START => $hsp->query->location->start,
				  Q_END   => $hsp->query->location->end,
		                  H_START => $hsp->hit->location->start,
				  H_END   => $hsp->hit->location->end,
				  SCORE   => $hsp->score, 
				  EVALUE  => $hsp->evalue,
				  H_SEQ   => $ss->seq(),
				  DIR     => $hsp->hit->location->strand,
				  });
  		  } 
  	      }
  	  }
      }
  }

  print Dumper($results);
