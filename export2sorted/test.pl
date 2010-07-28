#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;
use IO::File;

my $fh = IO::File->new();
$fh->open('Mash1_TC_Input/GA002-SR-R00027/s_7_filtered.txt');

for (my $i=1; $i<100; $i++){
    my @row = split "\t", <$fh>;
    die Dumper \@row  if $row[14] =~/\:/; 
}

$fh->close;
