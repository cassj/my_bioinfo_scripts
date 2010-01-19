#!/usr/bin/perl

use strict;
use warnings;
use Spreadsheet::ParseExcel;

use IO::File;

my $filename = shift @ARGV;
chomp $filename;

my $dir = shift @ARGV;
$dir = "." unless defined $dir;
chomp $dir;

die "File $filename not found" unless -e $filename;
die "File $filename not readable" unless -r $filename;

die "Dir not found" unless -d $dir;
die "Dir not writeable" unless -w $dir;

my $parser   = Spreadsheet::ParseExcel->new();
my $workbook = $parser->Parse($filename);

for my $worksheet ( $workbook->worksheets() ) {
 
  my $name = $worksheet->{Name};
  $name =~s/ /_/g;
  my $outfile = "$dir/$name.csv";
  
  my ( $row_min, $row_max ) = $worksheet->row_range();
  my ( $col_min, $col_max ) = $worksheet->col_range();

  next unless ($row_max>0 && $col_max>0);
  
  my $fh = new IO::File;
  $fh->open("> $outfile") or die "Can't open $outfile for writing: $!";

  for my $row ( $row_min .. $row_max ) {
    my @row;
    for my $col ( $col_min .. $col_max ) {
      my $cell = $worksheet->get_cell( $row, $col );
      my $value = $cell ? $cell->unformatted() : '';
      $value = '"'.$value.'"' if ($cell && ($cell->{Type} eq 'Text'));
      push @row, $value;
    }
    print $fh join ",", @row;
    print $fh "\n";
  }

  $fh->close;
}



