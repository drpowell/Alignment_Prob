#!/usr/bin/perl -w

use Markov_gen;

my $alphabet = [qw(a t g c)];

my $mm = new Markov_gen(1, $alphabet);
$mm->randModel();
$mm->model_power(10);

my $s1 = $mm->gen_sequence(10);
my $s2 = $mm->mutate($s1, 100);

print "s1=$s1\ns2=$s2\n";
print "H(s1) = ",$mm->entropy($s1),"\n";
print "H(s2) = ",$mm->entropy($s2),"\n";

print $mm->as_string(),"\n";
