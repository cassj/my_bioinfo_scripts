#!/usr/bin/perl

use strict;
use warnings;

use DBI;
use Data::Dumper;
use IO::File;

my @tables = qw/chr1   chr11  chr13  chr15  chr17  chr19  chr3   chr5   chr7   chr9   chrX chr10  chr12  chr14  chr16  chr18  chr2   chr4   chr6   chr8   chrMT  chrY /;


my $dbargs = {AutoCommit => 0,
	      PrintError => 1};

my $dbh = DBI->connect("dbi:SQLite:dbname=mmusculus_ncbi36_ensemblknowngene.db","","",$dbargs);

my $fh = new IO::File;
$fh->open("> mmusculus_ncbi36_ensemblknowngene.txt") or die "can't open file";

print $fh "RowNames\tChr\tGenomeStart\tGenomeEnd\tGenomeMid\tFeatureStart\tFeatureEnd\tStrand\tDescription\tEnsemblGeneID\n";

foreach (@tables){
    my $sth = $dbh->prepare("SELECT * FROM $_");
    $sth->execute();
    if ($dbh->err()) { die "$DBI::errstr\n"; }
    
    while ( my @row = $sth->fetchrow_array ) {
	my $line = join "\t",@row;
	print $fh "$line\n";
    }   

}

$fh->close;
$dbh->disconnect();

