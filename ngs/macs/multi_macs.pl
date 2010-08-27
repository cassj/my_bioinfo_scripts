#!/usr/bin/perl

use strict;
use warnings;

use IO::File;
use Cwd;

my $cmd = 'macs';
my $treatment = $ARGV[0];
my $control = $ARGV[1];
my $outdir = $ARGV[2] || '.';

my $genome_size = 2700000000; #this is the macs default anyway

#variables:
my @bandwidth = (100, 200, 300);
my @mfold = ('10,30';
my @pvalue = (0.001);

my $count = 1;

my $readme = new IO::File;
$readme->open(">$outdir/README") or die "Can't open file README for writing";

my $dir = getcwd;


#for each combinations of variable parameters, run macs.
foreach my $bw (@bandwidth){
    foreach my $mfold (@mfold){
        foreach my $pvalue (@pvalue){

            print "Running for BW: $bw ; MFOLD: $mfold ;  PVAL: $pvalue\n";
	    my $cmd = "macs -t '$dir/$treatment' --gsize $genome_size --bw $bw --mfold $mfold --pvalue $pvalue --format BED";
	    $cmd .= " -c '$dir/$control'" if $control;

            print $readme "$outdir/run$count\t$cmd\n";

            system("mkdir $outdir/run$count");
            chdir("$outdir/run$count");
            my $out = `$cmd 2>&1`;
            my $outfh = new IO::File;
            $outfh->open(">out.txt");
            print $outfh $out;
            $outfh->close;
            system("R --vanilla < NA_model.r 2>&1");
            chdir($dir);

            #pack up the results so we don't run out of space.
            #system("tar -czf $outdir/run$count.tgz $outdir/run$count");
            #system("rm -Rf $outdir/run$count");

            $count++;

        }
    }
}

$readme->close;



