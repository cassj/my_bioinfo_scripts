#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;

#############
#
# gtf2bed.pl
#
# Converts GTF files to BED files.
#
############

# should get this from cli.
print "track name=Transcript description='Converted from GTF' useScore=1\n";

while (my $line = <>) {
  chomp $line;
  my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attributes) = split "\t", $line;
  my %attributes = map {parse_att($_)} split ";", $attributes;

  my $name = $attributes{gene_id}.'_'.$attributes{transcript_id};
  $name =~ s/\"//g;
  print "$seqname\t$start\t$end\t$name\t$score\t$strand\n"

}

sub parse_att {
  my $att_string = shift;
  chomp $att_string;
  $att_string =~s/^\s+//;
  my ($key, $value) = split /\s/, $att_string;
  return ($key, $value);
}
