#!/usr/bin/perl -w
#
# Split a file with two sequences on the first two lines
# into two separate files in GENBANK format.  For use with
# the prss or ssearch programs

my($f) = shift;

my($str1,$str2);
open F, "< $f" || die "Can't read '$f'";
$str1 = <F>;
$str2 = <F>;

my $f1 = $f.".1";
my $f2 = $f.".2";

open(F,"> $f1") or die "Can't create tmp file";
print F ">GEN_IN1\n$str1\n";
close(F);
open(F,"> $f2") or die "Can't create tmp file";
print F ">GEN_IN2\n$str2\n";
close(F);
