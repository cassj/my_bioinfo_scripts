#!/usr/bin/perl

use strict;
use warnings;

use IO::File;

use Data::Dumper;

my $infile = new IO::File;
my $filename = $ARGV[0];
$infile->open("<$filename") or die "Can't open file $filename for reading";

# positions are 1-based (I checked with the UCSC browser, which is also 1-based despite 
# internal formats being 0-based).
# bed files are 0-based
# Positions are given on +ve strand regardless of tag strand
while (<$infile>){
    my ($tag_seq, $aln_score, undef, $chr_pos, $chr_dir, undef, undef) = split /\s/;
    my ($chr, $start) = split ':',  $chr_pos  ;
    my $length = length $tag_seq;
    my $end = $start+$length-1;
    $start = $start-1;
    $chr_dir = $chr_dir eq 'F' ? '+':'-';
    print "$chr\t$start\t$end\t$tag_seq\t$aln_score\t$chr_dir\n";
}

