#!/usr/bin/perl

use strict;
use warnings;

use IO::File;
use Cwd;

my $cmd = 'macs';
my $R = 'R';

my $treatment = $ARGV[0];
my $control = $ARGV[1];
chomp $treatment;
chomp $control;

#my $treatment = '../Mash1IP.bed';
#my $control = '../Mash1Input.bed';

my $genome_size = 2394590051; #70% of mm9, just a guess.

#variables:
my @bandwidth = (100, 200, 300);
my @mfold = (10,20,30);
my @pvalue = (0.001);

my $count = 1;

my $readme = new IO::File;
$readme->open('>README') or die "Can't open file README for writing";

#for each combinations of variable parameters, run macs.
foreach my $bw (@bandwidth){
    foreach my $mfold (@mfold){
        foreach my $pvalue (@pvalue){
            print "Running for BW: $bw ; MFOLD: $mfold ;  PVAL: $pvalue\n";

	    my $cmd = "macs -t '$treatment'";
	    $cmd .= " -c '$control'" if $control;
	    $cmd .= " --gsize $genome_size --bw $bw --mfold $mfold --pvalue $pvalue --format BED --wig";
            print $readme "run$count\t$cmd\n";
	    system("mkdir run$count");
            my $dir = getcwd;
            chdir("run$count");
            my $out = `$cmd  2>&1`;
            my $outfh = new IO::File;
            $outfh->open(">out.txt");
            print $outfh $out;
            $outfh->close;
            system("$R --vanilla < NA_model.r 2>&1");
            chdir($dir);

            #pack up the results so we don't run out of space.
            system("tar -czf run$count.tgz run$count");
            system("rm -Rf run$count");

            $count++;

        }
    }
}

$readme->close;
