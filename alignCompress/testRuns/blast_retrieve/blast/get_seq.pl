#!/usr/bin/perl -w

use Data::Dumper;

use Bio::DB::GenBank;
my $gb = new Bio::DB::GenBank();
#my $seq = $gb->get_Seq_by_id('AE001366.1');

my $seq;
my $id = shift;


if ($id =~ /^\d+$/) {
  $seq = $gb->get_Seq_by_gi($id);
} elsif ($id =~ /\./) {
  $seq = $gb->get_Seq_by_version($id);
} else {
  $seq = $gb->get_Seq_by_id($id);
}




print Dumper($seq);

#print "\n\nSEQUENCE\n";
#print $seq->seq()
