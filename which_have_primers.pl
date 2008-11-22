#!/usr/bin/perl

use strict;
use warnings;

use IO::File;

my $in_fh = new IO::File;
$in_fh->open("<agilent_chipchip_transcripts_mapping.csv") or die "can't open  agilent_chipchip_transcripts_mapping.csv for reading";

my $out_fh = new IO::File;
$out_fh->open(">coc_sites_with_primers.csv");

my $row = <$in_fh>;
chomp $row;
print $out_fh "$row\tHasPrimer\n";
while($row  = <$in_fh>){
  
  #ditch any quotes.
  chomp $row;
  $row =~ s/\"//g;
  my ($identifier, $trsc, $description,$TSS, $chr, $strand, $trans_start, $trans_end, $dist,$pos)=split("\t",$row);
  my $primer_file = "primers/Primer_Report_$trsc.txt";
  my $primer = -e $primer_file ? 'YES' : 'NO';
  print $out_fh "$row\t$primer\n";
}

$in_fh->close;
$out_fh->close;
