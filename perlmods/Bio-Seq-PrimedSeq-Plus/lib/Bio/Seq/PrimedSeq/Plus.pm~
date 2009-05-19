package Bio::Seq::PrimedSeq::Plus;

use strict;
use warnings;

our $VERSION = '0.01';

=head1 NAME 

  Bio::Seq::PrimedSeq::Plus - perl (bioperl) class to store a primed sequence along with information about BLAST results for the primers and thermodynamic parameters for both the primers and the amplicon

=head1 SYNOPSIS

  Extends Bio::Seq::PrimedSeq, which isa SeqFeature.

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


use base 'Bio::Seq::PrimedSeq';

#no, we'll use a local one
#use Bio::Tools::Run::RemoteBlast;
use Bio::Tools::Run::UNAFold::HybridSSMin;
    

=head2 new

  As for Bio::Seq::PrimedSeq

=head2 coerce

  A number of bioperl tools, notable Primer3, return 
  PrimedSeq objects. This method allows you to coerce 
  them into PrimedSeq::Plus objects
 

=cut

sub coerce{
  my ($class, $obj) = @_;

  #check we've got something approaching valid
  die "No object to coerce" unless $obj;
  die "Wrong type. ".ref($obj)." should be Bio::Seq::PrimedSeq objects."
    unless  (ref($obj) eq 'Bio::Seq::PrimedSeq'); 
  
  #rebless it
  return bless($obj, $class);

}


=head2 mfe_left_primer

  If called with no params, returns the results of the last
  call to the Unafold MFE folding algorithm hybrid-ss-min

  I called with params, will attempt to run hybrid-ss-min
  on the left primer sequence, passing any parameters to the 
  mfe algorithm.

=cut

sub mfe_left_primer{

  my $self = shift;
  my %params = @_;
  
  #if we have any arguments at all, redo mfe calcs
  if (%params){

    #save the params
    $self->left_mfe_params(\%params);
    
    #grab an hsmin processor and pass it params
    my $hsmin = Bio::Tools::Run::UNAFold::HybridSSMin->new(%params);
    $hsmin->verbose($self->verbose);       
    
    my $left_p = $self->get_primer('-left_primer');
    
    #Primer3 gives back IDs with spaces in, which causes
    #hybrid-ss-min to go boom.
    my $new_name = $left_p->display_id;
    $new_name =~ s/\s/_/g;
    $left_p->display_id($new_name);
    $hsmin->seq_obj($left_p->seq);
    $self->{mfe_left} = $hsmin->run;

  }

  return $self->{mfe_left};
}


=head2 left_mfe_params

  Returns a hashref of the params used for the last
  run of MFE. Or undef if it hasn't yet been run.

=cut

sub left_mfe_params{

  my ($self,$params) = @_;
  
  if (defined $params ){
    $self->{left_mfe_params} = $params;
  }
  return $self->{left_mfe_params};
	    
}


=head2 mfe_right_primer

  As for mfe_left_primer, but on the right primer

=cut

sub mfe_right_primer{
  my $self = shift;
  my %params = @_;
  
  #if we have any arguments at all, redo mfe calcs
  if (%params){
    
    $self->right_mfe_params(\%params);
    
    #grab an hsmin processor and pass it params
    my $hsmin = Bio::Tools::Run::UNAFold::HybridSSMin->new(%params);
    $hsmin->verbose($self->verbose);       
    
    my $right_p = $self->get_primer('-right_primer');
    
    #Primer3 gives back IDs with spaces in, which causes
    #hybrid-ss-min to go boom.
    my $new_name = $right_p->display_id;
    $new_name =~ s/\s/_/g;
    $right_p->display_id($new_name);
    $hsmin->seq_obj($right_p->seq);
    $self->{mfe_right} = $hsmin->run;

  }

  return $self->{mfe_right};

}



=head2 right_mfe_params

  Returns a hashref of the params used for the last
  run of MFE. Or undef if it hasn't yet been run.

=cut

sub right_mfe_params{
  my ($self,$params) = @_;
  if (defined($params)){
    $self->{right_mfe_params} = $params;
  }
  return $self->{right_mfe_params};
	    
}


=head2 mfe_amplicon

 As for mfe_left_primer, but for the entire amplicon.

=cut

sub mfe_amplicon{
  my $self = shift;
  my %params = @_;
  
  #if we have any arguments at all, redo mfe calcs
  if (%params){
    
    $self->amp_mfe_params(\%params);


    #grab an hsmin processor and pass it params
    my $hsmin = Bio::Tools::Run::UNAFold::HybridSSMin->new(%params);
    $hsmin->verbose($self->verbose);       
    
    my $amp = $self->amplicon;
    
    #Primer3 gives back IDs with spaces in, which causes
    #hybrid-ss-min to go boom.
    my $new_name = $amp->display_id;
    $new_name =~ s/\s/_/g;
    $amp->display_id($new_name);
    $hsmin->seq_obj($amp);
    $self->{mfe_amp} = $hsmin->run;
  }

  return $self->{mfe_amp};
}



=head2 amp_mfe_params

  Returns a hashref of the params used for the last
  run of MFE. Or undef if it hasn't yet been run.

=cut

sub amp_mfe_params{
  my ($self,$params) = @_;
  if (defined($params)){
    $self->{amp_mfe_params} = $params;
  }
  return $self->{amp_mfe_params};
	    
}

=head2 mfe_tms_under

  if($obj->mfe_tms_under(threshold=>60, folding_temp=>37)){...}

  Checks to see if the Tm values given by the MFE folding of
  both primers and the amplicon are all over the given threshold 
  at the given temperature.

  Warns and returns undef if any of the Tms are undefined 
  (ie if the folding hasn't been done yet). 
 
  Otherwise returns boolean 1 or 0.

  if no temp value is passed as an argument, the tmin
  from the mfe call is used.

=cut

sub mfe_tms_under{
  my ($self, %params) = @_;
  
  $self->throw( "No threshold given" ) unless $params{threshold};

  unless ($self->mfe_left_primer){
    warn "MFE not yet run on left primer";
    return undef;
  }
  unless ($self->mfe_right_primer){
    warn "MFE not yet run on right primer";
    return undef;
  }
  unless ($self->mfe_amplicon){
    warn "MFE not yet run on amplicon";
    return undef;
  }
 
  #if a temp is specified, check it is valid.
  if (defined $params{temp}){
    die "left primer MFE not run at that temperature" 
      unless defined $self->mfe_left_primer->tm->{$params{temp}};
    die "right primer MFE not run at that temperature" 
      unless defined $self->mfe_right_primer->tm->{$params{temp}};
    die "amplicon MFE not run at that temperature" 
      unless defined $self->mfe_amplicon->tm->{$params{temp}};
  }

  my $l_temp = $params{temp} || $self->mfe_left_primer->args->{tMin};
  my $r_temp = $params{temp} || $self->mfe_right_primer->args->{tMin};
  my $a_temp = $params{temp} || $self->mfe_amplicon->args->{tMin};
  
  my $l_tm =$self->mfe_left_primer->Tm->{$l_temp};
  my $r_tm =$self->mfe_right_primer->Tm->{$r_temp};
  my $a_tm =$self->mfe_amplicon->Tm->{$a_temp};
  
  
  return 1 if ($l_tm < $params{threshold} 
	       and $r_tm < $params{threshold}
	       and $a_tm < $params{threshold}
	      );
  
  return 0;  
}



=head2 _blast_defaults

  An internal method that returns a hash of blast params
  which have sensible defaults. (probably should be based on 
  type of sequence, primer design task and so on.)

=cut

sub _blast_defaults{
  my $self = shift;
  return  (
	   '-prog' => 'blastn',
	   '-data' => 'refseq_genomic',
	   '-readmethod' => 'SearchIO',		   
	   '-expect'=>10,
	  );
}


=head2 blast_left_primer

  If called with no arguments, returns the results of the last
  blast on the left primer. This may be undef if a blast has 
  yet to be run.

  If called with run=>true, it will attempt to run a remote 
  blast on the left primer, passing any other arguments to the 
  remote blast interface (see the documentation for 
  Bio::Tools::Run::RemoteBlast

  There are sensible (ymmv) defaults, so you can just do
   $obj->blast_left_primer(run=>1)

=cut

sub blast_left_primer{

  my ($self,%params) = @_;
  
  if( $params{run}){
    die "Trying to blast with no species defined" unless $self->species;
    
    delete $params{run};
    
    #set some sensible defaults, o/r with @_.
    #can't set readmethod cos other bits of the class expect it.
    %params = (
	       $self->_blast_defaults,
	       %params,
	      );

    my $factory = Bio::Tools::Run::RemoteBlast->new(%params);
    $Bio::Tools::Run::RemoteBlast::HEADER{'ENTREZ_QUERY'} = $self->species.'[orgn]';
   $Bio::Tools::Run::RemoteBlast::HEADER{'WORD_SIZE'} = '7';

#think this happens automatically when length < 30.
#    $Bio::Tools::Run::RemoteBlast::HEADER{'MEGABLAST'} = 'on';

    #Blast a sequence against a database:
    my $r = $factory->submit_blast($self->{left_primer}->seq);
    
    #grab result ids and retrieve the results
    foreach my $rid ( $factory->each_rid ) {
	
      # save the id
      $self->left_blast_rid($rid);
      
      my $rc;
      
      # loop until you get a response.
      warn "Blasting left primer ";
      while (!ref($rc)){
	  $rc = $factory->retrieve_blast($rid);
	  print '.' if $rc==0;
	  #need to flush stdout

	  #yell up if you get an error.
	  warn "BLAST received an error with rid: $rid, I'll keep trying, but if this keeps happening you may want to stop the process with ctrl-c " if $rc < 0;
	  #otherwise wait and try again.
	  sleep 10;
	}
	print "\n";
	

      # only one seq, one iteration so:
      $self->{blast_left_primer} = $rc->next_result;
      
    }
  }
  
  return $self->{blast_left_primer};
  
}


=head2 left_blast_rid
 
 returns the Retrieval ID of the last blast job on the
 left primer

=cut

sub left_blast_rid {
  my $self = shift;
  $self->{left_blast_rid} = shift if @_;
  return $self->{left_blast_rid};
}



=head2 left_blast_params

  A get method for the set of parameters passed to the
  last blast to be run on the left primer.

  Returns results as a list of key, value pairs which can 
  be put into a hash or handed direct to another blast call.

  my %params = $obj->left_blast_params
  $obj->blast_right_primer($obj->left_blast_params);


=cut 

sub left_blast_params{

 my $self = shift;

 if (@_){
   my $ps = {@_};
   $self->{left_blast_params} = $ps;
 }

 return @{$self->{left_blast_params}};

}

=head2 blast_right_primer

  As for blast_left_primer, but for the right primer.

=cut

sub blast_right_primer{

  my ($self,%params) = @_;
  
  if( $params{run}){
    die "Trying to blast with no species defined" unless $self->species;
    
    delete $params{run};
    
    #set some sensible defaults, o/r with @_.
    #can't set readmethod cos other bits of the class expect it.
    %params = (
	       $self->_blast_defaults,
	       %params,
	      );

    my $factory = Bio::Tools::Run::RemoteBlast->new(%params);
    $Bio::Tools::Run::RemoteBlast::HEADER{'ENTREZ_QUERY'} = $self->species.'[orgn]';
    $Bio::Tools::Run::RemoteBlast::HEADER{'WORD_SIZE'} = '7';
    

    #Blast a sequence against a database:
    my $r = $factory->submit_blast($self->{right_primer}->seq);
    
    #grab result ids and retrieve the results
    foreach my $rid ( $factory->each_rid ) {
	
      # save the id
      $self->right_blast_rid($rid);
      
      my $rc;
      
      # loop until you get a response.
      warn "Blasting right primer " ;
      while (!ref($rc)){
	  $rc = $factory->retrieve_blast($rid);
	  print '.' if $rc==0;
	  #yell if you get an error.
	  warn "BLAST received an error with rid: $rid, I'll keep trying, but if this keeps happening you may want to stop the process with ctrl-c " if $rc < 0;
	  #otherwise wait and try again.
	  sleep 10;
	}
	print "\n";
	

      # only one seq, one iteration so:
      $self->{blast_right_primer} = $rc->next_result;
      
    }



  }
  
  return $self->{blast_right_primer};
  
}



=head2 right_blast_rid
 
 returns the Retrieval ID of the last blast job on the
 right primer

=cut

sub right_blast_rid{
  my $self = shift;
  $self->{right_blast_rid} = shift if @_;
  return $self->{right_blast_rid};
}



=head2 right_blast_params

  As for left_blast_params, but for the last blast run on
  the right primer.

=cut


sub right_blast_params{

 my $self = shift;

 if (@_){
   my $ps = {@_};
   $self->{right_blast_params} = $ps;
 }

 return @{$self->{right_blast_params}};

}



=head2 blast_amplicon

  As for blast_left_primer, but for the whole amplicon
  sequence (including primers)

  Probably not really necessary, as blasting the primers
  should be enough to pick up any problem sequences but
  included for completeness.

=cut

sub blast_amplicon{

}


=head2 amplicon_blast_rid
 
 returns the Retrieval ID of the last blast job on the
 amplicon

=cut

sub amplicon_blast_rid{
  my $self = shift;
  $self->{amp_blast_rid} = shift if @_;
  return $self->{amplicon_blast_rid};
}




=head2 amplicon_blast_params

  As for left_blast_params, but for the last blast run
  on the amplicon.

=cut

sub amplicon_blast_params{

 my $self = shift;

 if (@_){
   my $ps = {@_};
   $self->{amp_blast_params} = $ps;
 }

 return @{$self->{amp_blast_params}};

}

=head2 blast_url

  Given a retrieval ID, will return the URL to retrieve
  that data using the NCBI web interface

  Makes no attempt to check whether the result is till 
  stored on the NCBI server.

=cut

sub blast_url{
  my ($self,$rid) = @_;
  die "can't make blast url without a retrieval ID." unless $rid;


  return 'http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?RID='.$rid.'&CMD=Web&PAGE_TYPE=BlastFormatting&SHOW_OVERVIEW=on&SHOW_LINKOUT=on&GET_SEQUENCE=yes&MASK_CHAR=2&MASK_COLOR=1&DESCRIPTIONS=100&ALIGNMENTS=100&NUM_OVERVIEW=100&GET_RID_INFO=on';

}


=head2 left_blast_url

  Returns the URL of the results page for the  last BLAST
  on the left primer. Or undef if there isn't one yet.

=head2 right_blast_url

  As above but for right primer

=head2 amplicon_url

  As above but for the amplicon

=cut

sub left_blast_url{
  my $self = shift;
  return undef unless $self->left_blast_rid;
  return $self->blast_url($self->left_blast_rid);
}

sub right_blast_url{
  my $self = shift;
  return undef unless $self->right_blast_rid;
  return $self->blast_url($self->right_blast_rid);
}

sub amplicon_blast_url{
  my $self = shift;
  return undef unless $self->amplicon_blast_rid;
  return $self->blast-url($self->amplicon_blast_rid);
}

=head2 blast_conflicts

  The blast results will tell you if a primer is hitting 
  a non-target region. A moderate level of homology with
  a non-target region may cause some competition for the 
  primer and it may not be an ideal choice, however if you 
  are constrained in your primer choice (for example if you
  are designing tiled primers) you may not wish to throw
  out a primer on this basis alone.

  Of far more serious concern is the case in which both
  primers have non-target hits close together, which could
  potentially result in amplification of a non-target region

  blast_conflicts will check the results of the left and right
  primer blasts and will return 1 if there are any conflicts, 0
  otherwise


  By default, a conflict is any instance where the primers hit
  closer than 4k bases. This can be changed with the min_separation 
  parameter although more than this is unlikely to cause problems
  in the time available for elongation in a normal pcr.

  my @conflicts = $obj->blast_conflicts(min_separation => '4000');

=cut

sub blast_conflicts{
  my ($self,%params) = @_;

  warn "testing for blast conflicts";
  
  die "You need to run blast_left_primer first" unless $self->blast_left_primer;
  die "you need to run blast_right_primer first" unless $self->blast_right_primer;
  
  my $min_sep = $params{min_separation} || 100000;

  #run blasts (or get cached results)
  my $left = $self->blast_left_primer;
  my $right = $self->blast_right_primer;

  #retrieve all hits
  my (@lefts, @rights);
  while (my $hit = $left->next_hit){push @lefts,$hit;}
  while (my $hit = $right->next_hit){push @rights,$hit;}

  #see if any of them conflict
  foreach my $left_hit (@lefts){
    foreach my $right_hit (@rights){
      if ($left_hit->strand == $right_hit->strand
	  && abs($right_hit->start - $left_hit->start) < $min_sep
	 ){

	#is it where we were expecting?
	warn "Left hit start ".$left_hit->start;
	warn "Right hit start ".$right_hit->start;
	warn "Official start ". $self->{left_primer}->start;
	return 1;
      } 
    }

    return 0;
  }
  return undef;
}



=head2 _genomic_dbs

  returns a hashref to a  map of species name to genomic db id.

=cut

sub _genomic_dbs{

    #'Homo sapiens [ORGN]';
#  return {
#	  mouse   => 'Mus musculus [ORGN]',
#	  rat     => 'Rattus norvegicus [ORGN]',
#	  chicken => 'Gallus gallus [ORGN]' ,
#	  macaca  => 'Macaca mulatta [ORGN]',
#	  chimp   => '',
#	  human   => 'Homo sapiens [ORGN]',
#	  dog     =>,
#	  cow     =>
#	 };
  
  return {
	  mouse     => '10090_genomic',
	  rat       => '10116_genomic',
	  chicken   => '9031_genomic',
	  macaca    => '9544_genomic',
	  chimp     => '9598_genomic',
	  human     => '9606_genomic',
	  dog       => '9615_genomic',
	  cow       => '9913_genomic',
	 };
}



=head2 repeat_dbs

 returns a hashref to a map of species name to repeat db id.

=cut

sub repeat_dbs{
  return {
	  mouse     => 'repeat_10090',
	  rat       => 'repeat_10116',
	  chicken   => 'repeat_9031',
	  macaca    => 'repeat_9544',
	  chimp     => 'repeat_9598',
	  human     => 'repeat_9606',
	  dog       => 'repeat_9615',
	  cow       => 'repeat_9913',
	 };
}



=head2 species

  Accessor for the species we are dealing with.
  Options are currently: human, mouse, rat, dog, cat
  cow, chimp, macaca

=cut

sub species{

  my $self = shift;

  my %valid_species;  
  foreach (qw (human mouse rat dog cat cow chimp macaca) )
    {$valid_species{$_}++}

  if (@_){
    my $sp = shift;
    if (defined $sp){
      $self->throw("Not a valid species: $sp") unless $valid_species{$sp};
    }
    $self->{species} = $sp;
  }

  return $self->{species};

}


1;
