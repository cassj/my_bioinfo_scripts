#!/usr/bin/perl

use strict;
use warnings;

use IO::File;

my $in_fh = new IO::File;
my $out_fh = new IO::File;

$in_fh->open('<Mash1CoC_Long_081104_raw.csv') or die "can't open infile";
$out_fh->open('>mm8_probes.csv') or die "can't open outfile";

#header
<$in_fh>;
#parse file
while (<$in_fh>){
    my (undef,undef,undef,undef,undef,$mm8_pos) = split(",",$_); 
    my ($chr, $start, $end) = $mm8_pos =~/^\"(.+):(\d+)-(\d+)\"$/;
    #if the probe is -ve strand, the positions are given in
    #-ve strand order, which is unconventional and causes
    #liftOver to get its knickers in a knot.
    #swap them over.
    ($start,$end) = ($end,$start) if $start > $end;
    print $out_fh "$chr $start $end\n";
}

$in_fh->close;
$out_fh->close;


#liftOver oldFile map.chain newFile unMapped

system('liftOver mm8_probes.csv /usr/local/share/liftoverfiles/mm8ToMm9.over.chain mm9_probes.csv unmapped');



#unmapped is empty, everything maps ok.





$in_fh->open('<Mash1CoC_Long_081104_raw.csv') or die "can't open infile";
my $in_fh2 = new IO::File;
$in_fh2->open("<mm9_probes.csv");	      
$out_fh->open('>Mash1CoC_Long_081104_raw_mm9.csv') or die "can't open outfile";

#header
my $row = <$in_fh>;
chomp($row);
print $out_fh "$row,mm9Position\n";

#parse file
while ($row = <$in_fh>){
    chomp($row);
    my $mm9 = <$in_fh2>;
    my ($chr,$start,$end) = split('\W',$mm9);
    print $out_fh "$row,$chr:$start-$end\n";
}

$in_fh->close;
$in_fh2->close;
$out_fh->close;
