#! /usr/bin/perl

use strict;
use warnings;
my $n, my $m;

my $infilename = shift(@ARGV);
my $outfilename = shift(@ARGV);

open(FILE, '<', $infilename) or die"$!";
open(OUT, '>', $outfilename) or die "$!";

my $firstline = <FILE>;
chomp($firstline);
if($firstline =~ /^\s*(\d+)\s(\d+)/) {
    $n = $1;
    $m = $2;
}
else {
    die "Not a valid metis file";
}

printf(OUT "%d\t%d\t%d\n", $n, $n, 2*$m); 


for(my $i=1; $i <= $n; $i = $i + 1) {
   my $line = <FILE>;
   chomp($line);
   my @list = split(/\s+/, $line);
   foreach my $j (@list) {
       if($j =~ /\S/){
	   printf(OUT "%s\t%s\t1\n", $i-1, $j-1);
       }
   }
}
