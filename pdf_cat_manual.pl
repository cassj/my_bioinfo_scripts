#!/usr/bin/perl

### pdf_pod.pl 
#
# Quick script to generate PDF files from Perl Module
# POD files (on Linux, not tested on Win)
#
###

use strict;
use warnings;

use File::Temp qw/tempfile tempdir/;
use File::Copy "cp";
use File::Find;
use File::Spec;
use PDF::API2;
use Paper::Specs;
use Getopt::Lucid qw( :all );

my  $opt = Getopt::Lucid->getopt( [
				   Switch("force_new|f"),
				   Param("pagesize|p")->default("a4"),
				   Switch("verbose|v"),
				  ] );


#Find your installed Catalyst Manual files

my @dirs;
my $mandir =  File::Spec->catfile('Manual.pm');
find (
      sub {
	print "$_\n";
	if(/.*Manual*/){
	  print "Helloooooo\n";
	  push @dirs, $_;
	}
      }, @INC
     );


use Data::Dumper;
die Dumper \@dirs;
#my @dirs = grep {/^.*$pattern$/}  @{$flh->locate($pattern)};


my $dir;

die "Catalyst::Manual not found. Are you sure it is installed?" if scalar(@dirs)==0;
if (scalar(@dirs)>1){
  print "Multiple installations of Catalyst::Manual found:\n";
  my $count=1;
  foreach (@dirs){
    print "[ $count ]\t $_\n";
    $count++;
  }
  print "Please enter the number of the location you wish to use:";
  my $ind = <>;
  print "\n";
  die "Invalid choice" unless $ind =~/^\d+$/ && $ind>0 && $ind<=scalar(@dirs); 
  $dir = $dirs[$ind-1];
}else{
  $dir = $dirs[0];
}

die "Can't read $dir" unless -r $dir;

unless (-w $dir){
  my $tmpdir = tempdir(UNLINK=>0);
  `cp -R $dir/* $tmpdir`;
  `chmod -R +w $tmpdir`;
  $dir = $tmpdir;
  warn "You don't have write permission to the Catalyst Manual directory, using $tmpdir instead";
}

##convert everything in Catalyst-Manual that is .pod to pdf
my @filenames;

find (
      sub {
	if (/\.pod$/ || /\.pm/){
	  warn "Converting $_ to pdf\n";
	  push @filenames, "$File::Find::dir/$_";
	  unless (-e "$File::Find::dir/$_.pdf" && !$opt->get_force_new){
	    system('pod2pdf '.
		    ' --output-file '.$_.'.pdf '.
		    ' --page-size '.$opt->get_pagesize.
		   ' '.$_
		  );
	  }
	}
      },
      ($dir)
     );


my $nfiles = scalar(@filenames);

if ($nfiles==0){
  die "\nNo .pod files to process. \n\nUse the -f switch to force recreation of pdfs from previously processed .pod files\n\n";
}

my @fileorder;
if ($nfiles == 1){
  $fileorder[0] = 1;
}
else{
  print "The following .pod files were found:\n";
  my $count = 1;
  foreach (@filenames){
    my ($filename) = /^$dir\/(.*)$/;
    print "[ $count ]\t$filename\n";
    $count++;
  }
    print "Please enter the numbers of the files you would like to include in the pdf, in the order in which you would like them, separated by spaces. For example: 1 4 2\n";
    
    my $fileorder = <>;
    chomp($fileorder);
    @fileorder = split '\s', $fileorder;
}
  
#check the file orders we've been given are valid.
foreach (@fileorder){
  die "Invalid file number $_" unless /^\d+$/ && $_>0 && $_< scalar(@filenames);
}

prFile("$dir/Manual.pdf");
prDoc("$filenames[$_-1].pdf") foreach @fileorder;

my @bookmarks;
my $pagenum=1;
foreach(@fileorder){
  my $file = $filenames[$_];
  my $this_pdf = PDF::API2->open("$file.pdf");
  my $section = $file =~ s/^.*\///;
  push @bookmarks, {text=> $section, act=> $pagenum};
  $pagenum += $this_pdf->pages;
  die $pagenum;
}

prBookmark( { text  => 'Document',
              close => 1,
              kids  => \@bookmarks } );


prEnd();

print "Done. Manual generated in $dir/Manual.pdf\n";
