package Bio::Annotation::RNAFold::HybridSSMin;


=head1 NAME 

  Bio::Annotation::RNAFold::HybridSSMin - Bioperl annotation object representing the results of HybridSSMin run

=head1 SYNOPSIS

 WARNING: This is the one to throw away. I'll get around to the rewrite shortly.

  Probably most of this stuff should be in the parent Bio::Annotation::RNAFold class. For now its in here cos I need it.

    Actually, maybe this shouldn't even be in the annotation namespace

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other Bioperl modules. Send your comments and suggestions preferably to one of the Bioperl mailing lists. Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General Discussion
  http://www.bioperl.org/MailList.html   - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track of the bugs and their resolution. Bug reports can be submitted via the web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR 

Cass Johnston

caroline.johnston@iop.kcl.ac.uk

=head1 CONTRIBUTORS

None yet but any offers of help would be gratefully accepted.

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceeded by a _

=cut

use strict;
use warnings;

our $VERSION = '0.01';

use Moose;
use Moose::Util::TypeConstraints;
use Bio::SeqIO;

#see http://search.cpan.org/~groditi/Moose-0.29/lib/Moose/Cookbook/FAQ.pod
#for extending non-moose classes and not screwing up constructors.
extends qw( 
	    Bio::Root::Root
	    Bio::AnnotationI
	    Moose::Object
	  );


subtype 'ReadableFile'
  => as 'Str'
  => where {-r $_}
  =>message {"Should be a valid, readable filename"};

subtype 'HashRefOfReadableFiles' 
  => as 'HashRef' 
  => where {
  my $test=1; 
  foreach my $file (values %$_)
    {$test=0 unless -r $file;} 
  return 1; } 
  => message{"All files should exist and be readable"};


=head2 args
 
 Returns a hash ref of the command line arguments
 passed to hybrid-ss-min to generate this result.
 
 Read Only (defined at object creation)

=cut

has args => (is => 'ro', isa => 'HashRef');


=head2 temps

 Returns an array ref of the temperatures at which 
 hybrid-ss-min was run. 

 Read Only (defined at object creation)

=cut

has temps => (is => 'ro', isa => 'ArrayRef');

=head2 prefix 

  Full path to the temp fasta file input to hybrid-ss-min
  which also serves as the prefix for the filename of 
  all the result files. 

  Although the file must exist and be readable, no effort
  is currently made to monitor the continued existance of
  the file.

=cut

has prefix => (is => 'ro', isa => 'ReadableFile');


=head2 dG

=cut

#need to work out how to limit these to being keyed by 
#$self->temps

has dH  => (is => 'ro', isa => 'HashRef');
has Tm  => (is => 'ro', isa => 'HashRef');
has dG  => (is => 'ro', isa => 'HashRef');
has dS  => (is => 'ro', isa => 'HashRef');
has img => (is => 'ro', isa => 'HashRefOfReadableFiles');




sub new {

  my $class = shift;

  my $obj = $class->SUPER::new(@_); 
  
  #not entirely sure what this is doing. 
  #makes Moose work. You're not really supposed to 
  #write your own constructors in Moose, but it's the 
  #only way to inherit from non-moose classes.
  my $self =  $class->meta->new_object(__INSTANCE__=> $obj, @_);

  return $self;
}


=head2 hash_tre

 This is required by Bio::AnnotationI
 I'm not entirely sure I see the point.

=cut 

sub hash_tree {
  my $self = shift;
  return {
	  temps => $self->temps,
	  dH    => $self->dH,
	  Tm    => $self->Tm,
	  dG    => $self->dG,
	  dS    => $self->dS,
	  img   => $self->img,
	  prefix => $self->prefix,
	 };
}


=head2 as_text

 This is required by Bio::AnnotationI

 Not really sure what a sensible text representation
 would be, but I've settled on the melting temp and 
 gibbs free energy of the MFE structure for each of the 
 folding temps tested

=cut

sub as_text{
  my $self = shift;

  my $txt;
  foreach (@{$self->temps}){
    $txt = "At $_: dG=".$self->dG->{$_}." Tm=".$self->Tm->{$_}.' ';
  }
  
  return $txt;

}


=head2 tagname 

  Default tagname for objects of this class when being
  added to Annotation::Collections.

  "mfe"

=cut

sub tagname{
  return "mfe";
}

1;
