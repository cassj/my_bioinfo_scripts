#!/usr/bin/perl

use strict;
use warnings;

use IO::File;
use Cwd;


#gah. Getopt, when I have time.
my $cmd = 'macs';
my $treatment = $ARGV[0];
my $control = $ARGV[1];
my $outdir = $ARGV[2] || getcwd;

my $genome_size = 2700000000; #this is the macs default anyway

#variables:
my @bandwidth = (100, 200, 300);
my @mfold = (10,20,30);
my @pvalue = (0.001);

my $count = 1;

my $readme = new IO::File;
$readme->open('>README') or die "Can't open file README for writing";

if (-e $outdir){
    opendir DIR, $outdir;
    if(grep !/^\.+$/, readdir(DIR)){
	die "Directory $outdir exists and is not empty.";
    }

}else{
    `mkdir $outdir`;
}


my $dir = getcwd;

#for each combinations of variable parameters, run macs.
foreach my $bw (@bandwidth){
    foreach my $mfold (@mfold){
        foreach my $pvalue (@pvalue){

            print "Running for BW: $bw ; MFOLD: $mfold ;  PVAL: $pvalue\n";

            print $readme "$outdir/run$count\tmacs -t '$dir/$treatment' -c '$dir/$control' --gsize $genome_size --bw $bw --mfold $mfold --pvalue $pvalue --format BED\n";

            `mkdir $outdir/run$count`;
            chdir("$outdir/run$count");
            
            my $out = `macs -t '$dir/$treatment' -c '$dir/$control' --gsize $genome_size --bw $bw --mfold $mfold --pvalue $pvalue --format BED  2>&1`;
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


