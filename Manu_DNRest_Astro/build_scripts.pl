#!/usr/bin/perl

## For Rory's data

use strict;
use warnings;

use Template;

my @cell_lines = ('Astro');
my $xpn_cols = {
		ESC => [qw(ev ev ev dn dn dn)],
		NS5 => [qw(ev ev ev dn dn dn)],
		Astro => [qw(dn dn dn dn ev ev ev ev)]
	       };

my @files = @ARGV;
chomp @files;

my $tt = Template->new({
    INTERPOLATE  => 1,
}) || die "$Template::ERROR\n";

foreach my $file (@files){
  foreach my $cell_line (@cell_lines){
   my ($filename) = $file =~ /.*\/(.+)\.tt2/;
    $tt->process($file, 
		 {cell_line => $cell_line,
		  xpn_cols => $xpn_cols,
		 }, "scripts/$cell_line/$filename")
      || die $tt->error(), "\n";
  }
}


