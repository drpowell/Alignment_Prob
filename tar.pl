#!/usr/bin/perl -w

use Archive::Tar;
my $tar = Archive::Tar->new;

my $fname = shift;
my $prefix = shift;

for my $f (@ARGV) {
  if ($f =~ /^(.*)=(.*)$/) {
    $tar->add_files($1);
    $tar->rename($1,$2);
  } else {
    $tar->add_files($f);
  }
}

$tar->write($fname, 1, $prefix);


