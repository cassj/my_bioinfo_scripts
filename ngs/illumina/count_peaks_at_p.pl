#!/usr/bin/perl

use strict;
use warnings;

use JSON;
use IO::File;

#am ignoring the FDR cos the calculation is a bit flaky.

my $fh = IO::File->new();
my @runs = 1..9;
my %counts;

my $p = shift @ARGV;
warn "counting peaks for p value $p";
my $score = -10*log($p)/log(10);
warn "converted p value to score $score";

foreach my $run (@runs){
  next unless -e "run$run/NA_peaks.xls";
  $fh->open("<run$run/NA_peaks.xls");

  #chr  start   end     length  summit  tags    -10*log10(pvalue)       fold_enrichment FDR(%)
  #chr1 3113438 3114067 630     326     23      34.03   4.46    33.63
  $counts{$run} = 0;

  while(<$fh>){
    next if /^#/;
    next if /^chr\s/;
    my @line = split "\t", $_;
    $counts{$run}++ if $line[6]>=$score;
  }
  $fh->close;
}

my $json = new JSON;
print $json->pretty->encode(\%counts);
