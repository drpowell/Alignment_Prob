#!/usr/bin/perl -w

use Archive::Tar;
my $tar = Archive::Tar->new;

my $fname = shift;
my $prefix = shift;

$tar->add_files(@ARGV);

$tar->write($fname, 1, $prefix);


