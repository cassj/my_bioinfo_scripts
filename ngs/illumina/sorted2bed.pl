#!/usr/bin/perl

use strict;
use warnings;
use IO::File;

use Data::Dumper;

my $filename = shift @ARGV;
die "Usage\n\tperl sorted2bed.pl s_X_sorted.txt > s_X_sorted.bed\n" unless $filename;
chomp $filename;
die "File $filename not found" unless -e $filename;
die "File $filename not readable" unless -r $filename;

my $fh = new IO::File;
$fh->open("< $filename") or die "Can't open file $filename for reading: $!";

my $count = 1;
while(my $line = <$fh>){
   warn "Line $count\n" if $count%1000 == 0;
   $count++;
   my @line = split "\t", $line;
   my $chr = $line[10];
   $chr =~ s/(.+)\.fa/$1/;
   #Illumina is 1-based, BED is 0-based
   my $start = $line[12]-1;
   my $read = $line[8];
   my $end = $start + length $read;
   my $strand = $line[13] eq 'F' ? '+': '-';
   my $score = $line[15];
   my $bedline = "$chr\t$start\t$end\t$read\t$score\t$strand\n";
   print $bedline;
}
$fh->close;

warn "Done";
