#!/usr/bin/perl -w

use Data::Dumper;

use Bio::DB::GenBank;
my $gb = new Bio::DB::GenBank();
#my $seq = $gb->get_Seq_by_id('AE001366.1');
my $seq = $gb->get_Seq_by_id($ARGV[0]);



print Dumper($seq);

print $seq->seq()
