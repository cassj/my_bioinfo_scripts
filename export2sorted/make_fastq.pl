#!/usr/bin/perl

#converts an illumina export.txt file to fastq
#Note that this doesn't do any score conversion. 
#make sure you hand the appropriate cmdline switch to bowtie
#either --solexa-quals or --solexa1.3-quals

use strict;
use warnings;

use IO::File;
use Data::Dumper;

my $fh = new IO::File;
my $filename = $ARGV[0];
die "Usage:\n\tmake_fastq.pl export.txt > export.fastq" unless $filename;
 
$fh->open($filename) or die "can't open file $filename: $!";

my $i=1;
while(<$fh>){
  my @row = split "\t", $_;
  
  my $filter = pop @row;
  #deal with windows or nix line endings. 
  $filter =~ s/\r?\n$//;
  #failed filter?
  next if ($filter eq 'N');

  #id is Machine-runnum-lane-tile-x-y
  my $id = join '-', @row[0..5];
  my $seq = $row[8];
  my $qual = $row[9];
  print "\@$id\n$seq\n+\n$qual\n";
}

$fh->close;

