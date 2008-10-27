=head1 NAME

Bio::Annotation::ExternalLocation - Annotate sequences with an external location

=head1 SYNOPSIS

  my $annot = Bio::Annotaion:ExternalLocation->new(
                -start             => $start,
                -end               => $end,
                -strand            => -1,
                -taxon             => $taxon,
                -authority         => 'NCBI',
                -coord_sys_type    => 'chromosome'
                -coord_sys_version => 37,
                -coord_sys_id      => 'X',
                -tagname           => $tagname,
  );


=head1 DESCRIPTION

Store an external location for a sequence, typically a location in a 
genome build (NCBI, UCSC Golden Path etc).
Implements Bio::AnnotationI and Bio::RangeI

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@biperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Cass Johnston

Email:  caroline-dot-johnston-at-iop-dot-kcl-dot-ac-dot-uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::Annotation::ExternalLocation;

use strict;
use Carp;

use base qw(Bio::Range Bio::AnnotationI);

sub new {
  my ( $class, @args) = @_;

  #set RangeI and AnnotationI values
  my ($self) = $class->SUPER::new(@args);

  #set a default tagname
  $self->tagname('ExternalLocation');

  #and sort out our args
  my %args = @args;
  foreach (keys %args) {
    #deal with the bioperl dashes.
    my ($method) = /^-?(.+)$/;
    #and set the value
    $self->$method($args{$_});
  }


  return $self;

}

=head2 authority

  Title   : authority
  Usage   : $auth = $loc->authority($new_auth)
  Function: Accessor for the authority / institution that 
            defines the accession codes of a co-ord system
            or that provides a gene build
  Returns : [String] The authority name
  Args    : Optionally,a new value for authority

=cut

sub authority {
  my ($self,$value) = @_;
  if( defined $value) {
    $self->{'authority'} = $value;
  }
  return $self->{'authority'};
}


=head2 taxon

  Title   : taxon
  Usage   : $taxon = $loc->taxon($taxon)
  Function: Accessor for the species
  Returns : An object of class Bio::Taxon
  Args    : Optionally,a new value for taxon

=cut

sub taxon {
  my ($self,$value) = @_;
  if( defined $value) {
    warn 'Taxon is not a Bio::Taxon object. This may cause problems.' unless $value->isa('Bio::Taxon');
    $self->{'taxon'} = $value;
  }
  return $self->{'taxon'};
}




=head2 url

  Title   : url
  Usage   : $url = $loc->url($new_url)
  Function: Accessor for a url associated with this annotation
  Returns : [String] URL associated with this annotation, or undef
  Args    : [String] Optional new URL

=cut

sub url {
  my ($self, $value) = @_;
  if( defined $value) {
    $self->{'url'} = $value;
  }
  return $self->{'url'};
}


=head2 coord_sys_type

  Title   : coord_sys_type
  Usage   : $coord_sys_type = $loc->coord_sys_type()
  Function: Accessor for the coordinate system type
  Returns : [String] The coordinate system type
  Args    : [String] Optional new coordinate system 
            type

=cut

sub coord_sys_type {
  my ($self, $value) = @_;
  if( defined $value) {
    $self->{'coord_sys_type'} = $value;
  }
  return $self->{'coord_sys_type'};
}


=head2 coord_sys_id

  Title   : coord_sys_id
  Usage   : $coord_sys_id = $loc->coord_sys_id()
  Function: Accessor for the coordinate system id
  Returns : [String] The coordinate system id (eg. a chromosome name)
  Args    : [String] Optional new coordinate system id

=cut

sub coord_sys_id {
  my ($self, $value) = @_;
  if( defined $value) {
    $self->{'coord_sys_id'} = $value;
  }
  return $self->{'coord_sys_id'};
}


=head2 coord_sys_version

  Title   : coord_sys_version
  Usage   : $coord_sys_version = $loc->coord_sys_version($new_version)
  Function: Accessor for the coordinate system version
  Returns : [String] The coordinate system version
  Args    : [String] Optional new coordinate system version

=cut

sub coord_sys_version {
  my ($self, $value) = @_;
  if( defined $value) {
    $self->{'coord_sys_version'} = $value;
  }
  return $self->{'coord_sys_version'};
}





=head1 Bio::AnnotationI Implementation

=head2 start

  Title   : start
  Usage   : $start = $loc->start($new_start)
  Function: Accessor for the start of the location 
            By convention, positions are given relative to the
            positive strand, so start is always less than end
            regardless of the strand of interest.
  Returns : [Int] The start position
  Args    : [Int] Optionally a new start position

=cut

sub start{

  my ($self, $val) = @_;
  if (defined $val){
    $self->{start} = $val;
  }
  return $self->{start};
}



=head2 end

  Title   : end
  Usage   : $end = $loc->end($new_end)
  Function: Accessor for the end of the location 
            By convention, positions are given relative to the
            positive strand, so start is always less than end
            regardless of the strand of interest.
  Returns : [Int] The end position
  Args    : [Int] Optionally a new end position

=cut

sub end{

  my ($self, $val) = @_;
  if (defined $val){
    $self->{end} = $val;
  }
  return $self->{end};
}




=head2 strand

  Title   : strand
  Usage   : $strand = $loc->strand($new_strand)
  Function: Accessor for the strand of the location 
            1 = Positive Strand
            -1 = Negative Strand
  Returns : [1 or -1] The strand of the sequence
  Args    : [1 or -1] Optionally a new value for strand

=cut

sub strand{

  my ($self, $val) = @_;
  if (defined $val){
    $self->{strand} = $val;
  }
  return $self->{strand};
}



=head2 length

  Title   : length
  Usage   : $len = $loc->length()
  Function: calculates the length of the sequence from 
            start and end values
  Returns : [Int] The length of the sequence or undef if
            either start or end are undefined
  Args    : none.

=cut

sub length{

  my $self = shift;
  return undef unless defined $self->start;
  return undef unless defined $self->end;
  return ($self->end - $self->start) + 1;
}





=head1 Bio::AnnotationI Implementation

=head2 as_text

  Title   : as_text
  Usage   : print $annot->as_text
  Function: Returns a text description of the DB location
  Returns : String
  Args    : none

=cut

sub as_text {
  my ($self) = shift;

  my $auth = $self->authority || 'Source_Unknown';
  my $species = defined $self->taxon ? $self->taxon->scientific_name : 'Species_Unknown';
  my $start = $self->start || '??';
  my $end = $self->end || '??';
  my $strand = $self->strand || '??';
  my $length = $self->length || '??';
  my $coord_type = $self->coord_sys_type || 'Unknown Coord System Type';
  my $coord_version = $self->coord_sys_version || 'Unknown Coord System Version';
  my $coord_id = $self->coord_sys_id || 'Unknown Coord System ID';
  
  return "Auth: $auth; Species: $species; Version: $coord_version; Type: $coord_type; ID: $coord_id; $start to $end ($length bases on strand $strand).";
}



=head2 hash_tree

  Title   : hash_tree
  Usage   : my $hashref =  $annot->hash_tree
  Function: Returns the location in the form of a hash
  Returns : HashRef
  Args    : none

=cut
sub hash_tree{
  my $self = shift;

  my $species = defined $self->taxon ? $self->taxon->scientific_name : '';

  my $start = $self->start || '??';
  my $end = $self->end || '??';
  my $strand = $self->strand || '??';
  my $length = $self->length || '??';
  my $coord_type = $self->coord_sys_type || 'Unknown Coord System Type';
  my $coord_version = $self->coord_sys_version || 'Unknown Coord System Version';
  my $coord_id = $self->coord_sys_id || 'Unknown Coord System ID';
  
  return {
	  authority         => $self->authority,
	  species           => $species,
	  coord_sys_type    => $self->coord_sys_type,
	  coord_sys_version => $self->coord_sys_version,
	  coord_sys_id      => $self->coord_sys_id,
	  start             => $self->start,
	  end               => $self->end,
	  strand            => $self->strand,
	  length            => $self->length,
	 };
}


=head tagname

  Title   : tagname
  Usage   : $annot->tagname($newval)
  Function: Get/set the tagname (default is GenomePosition)
  Returns : String
  Args    : optional new value for tagname

=cut
sub tagname {
  my ($self, $val) = @_;
  if (defined $val){
    $self->{tagname} = $val;
  }
  return $self->{tagname};
  
}


1;
