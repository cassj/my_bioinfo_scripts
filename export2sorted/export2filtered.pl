#!/usr/bin/perl

#converts an illumina export.txt file to filtered.txt
#by removing lines which fail filtering or have no 
#chromosomal mapping and removing the final (filter Y/N)
#column. Also swaps windows newline for unix.

use strict;
use warnings;

use IO::File;
use Data::Dumper;

my $fh = new IO::File;
my $filename = $ARGV[0];
die "Usage:\n\texport2filtered.pl export.txt > filtered.txt" unless $filename;
 
$fh->open($filename) or die "can't open file $filename: $!";

while(<$fh>){
  my @row = split "\t", $_;
  my $filter = pop @row;
  #deal with windows or nix line endings. 
  $filter =~ s/\r?\n$//;
  #failed filter?
  next if ($filter eq 'N');
  #failed mapping?
  next unless ($row[13]);  
  print join "\t", @row;
  print "\n";
}

$fh->close;

