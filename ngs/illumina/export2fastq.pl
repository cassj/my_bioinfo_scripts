#!/usr/bin/perl
use strict;
use warnings;

while(<>){
  chomp;
  my @parts = split /\t/;
  $parts[8]=~s/\./N/g;
  print "@","$parts[0]:$parts[2]:$parts[3]:$parts[4]:$parts[5]#$parts[6]/$parts[7]\n";
  print "$parts[8]\n";
  print "+\n";
  print "$parts[9]\n";
}

