#!/usr/bin/perl;

use strict;
use warnings;
use Data::Dumper;

sub remove_suffix{
    my $filename = shift;
    if ($filename =~/(.+)\.tt2?/ ){
	rename $filename,  $1;
    }
}

sub do_dir{
    my $dir = shift;
    opendir DIR, $dir;
    my @files = grep{!/^\.+$/} readdir DIR; 

    foreach (@files){
	
	if (-d "$dir/$_"){
	    &do_dir("$dir/$_");
	}
	else{
	    &remove_suffix("$dir/$_");
	}
	
    }
}

&do_dir($ARGV[0]);
