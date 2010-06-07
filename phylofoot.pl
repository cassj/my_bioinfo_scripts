#!/usr/bin/perl

# phylofoot.pl
#   by Boris Lenhard
#   (variously faffed with by CassJ. Original from 
#   http://tfbs.genereg.net/)

use strict;
use warnings;

use Getopt::Long; 
use Pod::Usage;
use TFBS::DB::FlatFileDir;

my ($database_dir, $alignment_file, $help);
my @matrix_IDs   = ();
my $conservation = 70;
my $threshold    = 80;
my $window       = 50;

GetOptions('help'              => \$help,
	   'alignment_file=s'  => \$alignment_file,
	   'database=s'        => \$database_dir,
	   'matrix-id:s'       => \@matrix_IDs,
	   'conservation:f'    => \$conservation,
	   'threshold-score:f' => \$threshold,
	   'window-size:i'     => \$window
	   );

if($help)  {
    pod2usage(-exitstatus=>0, -verbose=>2);
}
elsif (!($alignment_file and $database_dir)) {
    pod2usage(1);
}

  # parse both the comma-separated lists and individually specified
  # matrix IDs:

@matrix_IDs = split (",", join(',',@matrix_IDs));


my $db = TFBS::DB::FlatFileDir->connect($database_dir);

my $matrixset;

unless (scalar @matrix_IDs) {
  $matrixset = $db->get_MatrixSet(-matrixtype=>"PWM");
}
else {
  $matrixset = $db->get_MatrixSet(-IDs        => \@matrix_IDs,
				  -matrixtype => "PWM");
}


# do the phylogenetic footprinting search of the sites

my $sitepairset = $matrixset->search_aln(-file=>$alignment_file,
					 -cutoff=>"$conservation",
					 -threshold=>"$threshold\%",
					 -windowsize => $window);

print $sitepairset->GFF;









__END__


=head1 NAME

phylofoot.pl - Phylogenetic footprinting example script

=head1 SYNOPSIS

./phylofoot.pl -a <alignment_file> -d <TFBS_matrix_dbase_dir> [other_options] 

=head1 OPTIONS

=over 8

=item B<-d   or --database>  <directory name>

REQUIRED: Name of the FlatFileDir database directory to 
use for retrieving matrices. 
A sample database directory examples/SAMPLE_FlatFileDir 
is available in TFBS distribution. 

=item B<-a   or --alignment-file> <alignment file name>

REQUIRED: Name of the pairwise alignment file in Clustal 
format.
A sample database directory examples/sample_alignment.aln 
is available in TFBS distribution. 

=item B<-m   or  --matrix-id> <list of matrix IDs>

OPTIONAL: ID of the matrix from the database to scan the 
alignment with.You can specify multiple matrices using 
multiple -m switches or a single comma-separated lists of 
IDs (NO spaces - e.g. -m M00001,M00021,N01921 ). You can 
use a script called examples/list_matrices.pl in TFBS distribution
to list information for all matrices in a matrix database 
of the FlatFileDir type.

DEFAULT: If no matrix IDs are specified, all matrices in 
the database are used for the search;

=item B<-w   or  --window-size> <integer value>

OPTIONAL: The width of sliding window for calculating
the conservation profile of the submitted pairwise alignment.

DEFAULT: If not specified, the default value is 50 (nucleotides).

=item B<-c   or  --conservation> <percent value>

OPTIONAL: Conservation cutoff (in percent) for a region of
multiple alignment to include detected conserved sites into output.

DEFAULT: If not specified, the default value is 70 (%).

=item B<-t   or  --threshold-score> <percent value>

OPTIONAL: Threshold score (in percent) for a matrix match to 
a subsequence.

DEFAULT: If not specified, the default value is 80 (%).

=back



=head1 DESCRIPTION

This is an example script that scans conserved regions of a 
pairwise DNA sequence alignment with a set of matrices form 
a flat file databases and produces GFF output. Its source code is 
meant to be studied by bioinformaticians who wish to learn how to 
use TFBS modules.


=cut

