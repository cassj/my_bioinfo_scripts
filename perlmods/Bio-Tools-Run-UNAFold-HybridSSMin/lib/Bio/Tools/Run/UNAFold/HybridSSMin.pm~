package Bio::Tools::Run::UNAFold::HybridSSMin;

=head1 NAME 

  Bio::Tools::Run::UNAFold::HybridSSMin - bioperl run wrapper for the UNAFold hybrid-ss-min program.

=head1 SYNOPSIS

  First attempt at a Bioperl Run wrapper for the 
  UNAFold hybrid-ss-min tool.

  Requires the installation of UNAFold from
  http://www.bioinfo.rpi.edu/applications/hybrid/download.php

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
use Tie::File;
use Cwd;
use Scalar::Util qw/weaken/;

use Bio::Annotation::RNAFold::HybridSSMin;

#see http://search.cpan.org/~groditi/Moose-0.29/lib/Moose/Cookbook/FAQ.pod
#for extending non-moose classes and not screwing up constructors.
extends qw( Bio::Root::Root 
            Bio::Tools::Run::WrapperBase
	    Moose::Object
	  );


use vars qw/@cmd_args $convert_bin $draw_bin/;

BEGIN{
  
  @cmd_args = qw/NA tmin tinc tmax sodium magnesium polymer suffix output prohibit force energyOnly noisolate mfold maxbp constraints basepairs circular allpairs maxloop nodangle simple prefilter/;
  
}

#really? surely these should be in the Vienna stuff?
#$convert_bin =  "/usr/local/share/ViennaRNA/bin/ct2b.pl";
#$draw_bin = "/usr/local/bin/RNAplot";  




=head1 hybrid-ss-min Parameters

get/set accessors are provided for the following parameters for hybrid-ss-min:

NA          : (RNA|DNA) Nucleic acid type. Default is RNA
tmin        : (int) Minimum temperature. Default 37
tinc        : (int) Temperature increment. Default 1
tmax        : (int) Maximum temperature. Default 37
sodium      : (float) Sodium ion concentration (molar). Default is 1
magnesium   : (float) Magnesium ion concentration (molar). Default is 0
polymer     : (boolean) Use salt corrections for polymers instead of oligos. Default is 1
suffix      : (string) Use energy rules with the suffix string. Overrides tmin, tinc, tmax, sodium, magnesium and polymer
output      : (string) prefix output files with output
prohibit    : (i,j,k) prohibit all basepairs in the helix from i,j to i+k-1, j-k+1. If j is 0, prohibit bases i to i+k-1 from pairing at all. If i is 0, prohibit bases j to j-k+1 from pairing at all. k defaults to 1 if not specified
force       : (i,j,k) force all basepairs in the helix from i,j to to i+k-1, j-k+1. If j is 0, forces bases i to i+k-1 to be double-stranded. If i is 0, forces bases j to j-k+1 to be double-stranded. k defaults to 1 if not specified
energyOnly  : (boolean) skips computation of structure and returns only dG
noisolate   : (boolean) prohibits all isolated base-pairs (helices of length 1).
mfold       : (P,W,MAX or boolean) causes hybrid-ss-min to perform multiple suboptimal tracebacks in the style of MFold. P indicates the percent sub-optimality to consider. only structures with energies within P% of the minimum will be output. W indicates the window size. A structure must have at least W base pairs that are each a distance of at least away from any basepair in a previous stucture. MAX represents an absolute limit on the number of structures computed. This option also tells hybrid-ss-min to generate p-num values and energy dot plot. 
maxbp       : (int) bases further apart than this value cannot form pairs. Default is no limit.
constraints : (filename) reads a list of constraints from a file. Constraints must be in the form "Pijk" or "Fijk". These are equivalent to specifying PROHIBIT or FORCE ijk. If not sepecifed, filename defaults to PREFIX.aux
basepairs   : (filename) reads a list of allowable helices from a file. Each helix consists of three whitespace-delimited numbers which specify the starting basepair and the length of the helix. When this option is used all basepairs except those in <filename> are prohibited from forming.
circular    : (boolean) treats sequences as circular rather than linear
allpairs    : (boolean) allows basepairs to form between any two nucleotides. Default is only Watson-Crick and Wobble pairings
maxloop     : (int) Maximum size of bulge/interior loops. Default is 30.
nodangle    : (boolean) removes single-base stacking from consideration
simple      : (boolean) makes the penalty for multibranch loops constant rather than affine
prefilter   : (val1 val2) Sets the prefilter to filter out all basepairs except those in groups of val2 adjacent base-pairs of which val1 can form. val2 is the same as val1 if unspecified. Default is 2 of 2.

Usage in all cases is:

 #get
 my $val = $obj->method;

 #set
 $obj->method($val);

 #set at object creation
 my $obj = Bio::Tools::Run::UNAFold::HybridSSMin->new(method1=>'value',method2=>'value2');

=cut

#hybrid-ss-min args
#readablefile is provided by Annotation::RNAfold::HybridSSMin. This is a bit of a mess really.

subtype 'PosInt' => as 'Int' =>  where { $_ > 0};
subtype 'PosOrZeroInt' => as 'Int' => where {$_ >= 0 };
subtype 'PosOrZeroFloat'=> as 'Num' => where {$_ >= 0};
subtype 'ijk' => as 'Str' => where {/^\d+,\d+,\d+$/} => message {"Should be in the form: i,j,k"};
subtype 'PWMax' => as 'Str' => where {/^\d+,\d*,?\d*$/ || /^[10]$/} => message {"Should be in the form P,W,Max"};
#subtype 'ReadableFile'=> as 'Str'=> where {-r $_}=>message {"Should be a valid, readable filename"};
subtype 'Prefilter' => as 'Str' =>where {/^\d+,\d+$/}=> message {"Should be in the form: val1,val2"};
subtype 'Dir'=> as 'Str' =>where {-d $_}=>message {"Should be a valid directory"};
subtype 'SeqObj' => as 'Object' => where {$_->isa('Bio::SeqI') && defined $_->seq && '' ne $_->seq }=> message {'Should be a Bio::SeqI implementation with a sequence defined'};
subtype 'NucAc' => as 'Str' => where {/^(D|R)NA$/i};

has 'NA' => (is => 'rw', isa => 'NucAc|Undef');
has 'tmin' =>(is => 'rw', isa   => 'PosOrZeroInt|Undef');
has 'tinc' =>(is => 'rw', isa => 'PosOrZeroInt|Undef');
has 'tmax' =>(is    => 'rw',isa   => 'PosOrZeroInt|Undef');
has 'sodium' => (is => 'rw', isa => 'PosOrZeroFloat|Undef');
has 'magnesium' => ( is => 'rw', isa => 'PosOrZeroFloat|Undef');
has 'polymer' => (is => 'rw', isa => 'Bool|Undef');
has 'suffix' => (is => 'rw', isa => 'Str|Undef');
has 'output' => (is => 'rw',isa => 'Str|Undef');
has 'prohibit' => (is => 'rw', isa => 'ijk|Undef');
has 'force' => (is => 'rw',isa => 'ijk|Undef');
has 'energyOnly' => (is => 'rw',isa => 'Bool|Undef');
has 'noisolate' => (is => 'rw', isa => 'Bool|Undef');
has 'mfold' => (is => 'rw', isa => 'PWMax|Undef');
has 'maxbp' => (is => 'rw',isa => 'PosInt|Undef');
has 'constraints' => (is => 'rw',isa => 'ReadableFile|Undef');
has 'basepairs' => (is => 'rw', isa=>'ReadableFile|Undef');
has 'circular' => (is => 'rw',isa => 'Bool|Undef');
has 'allpairs' => (is => 'rw', isa => 'Bool|Undef');
has 'maxloop' => (is => 'rw',isa => 'PosInt|Undef');
has 'nodangle' => (is => 'rw',  isa => 'Bool|Undef');
has 'simple' => (is => 'rw',isa => 'Bool|Undef');
has 'prefilter' => (is => 'rw',isa => 'Prefilter|Undef');


# other stuff

has 'program_name' => (is => 'rw',isa => 'Str', default => 'hybrid-ss-min');
has 'program_dir' => (is => 'rw',isa => 'Dir');
has 'seq_obj' => (is => 'rw', isa => 'SeqObj', predicate => 'has_seq' );
has 'plot_dir' => (is => 'rw', isa => 'Dir');




=head1 WrapperBase Methods

  These methods are directly inherited from WrapperBase and not overriden by this class. 
  Please see that module's perldoc for details. 

=head2 arguments

=head2 no_param_checks

=head2 save_tempfiles

=head2 outfile_name

=head2 tempdir

=head2 cleanup

=head2 io

=head2 executable


=head1 Overridden WrapperBase methods 

=head2 program_dir

  Usage       : my $dir = $obj->program_dir; $obj->program_dir($dir);
  Function    : get/set the directory in which the executable is stored.
  Default     : '', which should be fine if the executable is in your PATH
  Note        : Defined but not implemented in WrapperBase.

=head2 program_name

  Usage       : my $name = $obj->program_name; $obj->program_name($name);
  Function    : get/set the name of the executable.
  Default     : 'hybrid-ss-min'
  Note        : Defined by not implemented in WrapperBase

=head1 CONSTRUCTOR

=head2 new

 Title     : new()
 Usage     : my $hybssmin = Bio::Tools::Run::UNAFold::HybridSSMin->new(seq=>$aseqobj);
 Function  :
 Returns   : An object of class Bio::Tools::Run::UNAFold::HybridSSMin
 Args      : All arguments are optional as they can be added later

=cut

sub new {

  my $class = shift;

  #call Bio::Root::Root constructor
  my $obj = $class->SUPER::new(@_); 
  
  #kludge moose and non-moose parent classes
  my $self =  $class->meta->new_object(__INSTANCE__=> $obj, @_);

  return $self;
}



=head2 run

 Usage    : $obj->run;
 Function : Runs hybrid-ss-min on the data in the object.
 Returns  : For now, a hash of results
 Args     : None
 Notes    : Most results will be returned in the resulting hash.
            Postscript images of the MFE structures will only be 
            returned if $obj->save_tempfiles is true. Once you've 
            done whatever you want to do with the plotfiles (in 
            $res->{images}->{temp}), for example saving them to 
            somewhere else, then you can remove the tempfiles by 
            calling $obj->cleanup. 

=cut


sub run{
  my $self = shift;

  #good to go?
  if ($self->_check_setup){
 
    # use tempfile to get a random filename.
    # keep it. cleanup will deal with it later.
    my ($tfh,$tempfile) = $self->io->tempfile(DIR=> $self->tempdir, UNLINK=>0);
    
    # use tempfile name for your output prefix.
    $self->output($tempfile);
    
    # now write your sequence to that file in fasta
    my $infile_io = Bio::SeqIO->new('-fh' => $tfh,
				    '-format' => 'fasta',
				   );
    $infile_io->write_seq($self->seq_obj);
    close $tfh;

    print "written file as $tempfile\n" if $self->verbose;

    # concatenate bits of cmd.
    my $cmd = $self->executable.' '.$self->arguments.' '.$tempfile;
    print "\nRunning: \n\t$cmd\n\n" if $self->verbose;
    
    # hybrid-ss-min outputs a load of crap. Redirect it:
    open(OLDOUT, ">&STDOUT"); #copy file descriptor
    open(STDOUT, ">> $tempfile.log") 
      or die "Can't redirect STDOUT to $tempfile.log";
    
    #run the cmd.
    my $status = system("$cmd");
    $self->throw($self->program_name." crashed: $? $cmd\n")
      unless ($status==0) ;

    close(STDOUT)                       
      or die "Can't close STDOUT: $!";    #close redirected fh
    open(STDOUT, ">&OLDOUT")            
      or die "Can't restore stdout: $!";  #reopen old one
    close(OLDOUT)                       
      or die "Can't close OLDOUT: $!";    #close copied descriptor

    #parse the results
    my $res = $self->_parse_results($tempfile);
  
    # for reasons I'm not sure of this isn't getting
    # called, My guess would be that DESTROY doesn't
    # play nicely with Moose, at least with inheritance
    $self->cleanup unless $self->save_tempfiles;

    return $res;

  }

}

##### TEMP SOLUTION : THIS NEEDS COMPLETELY RE-DONE

sub _parse_results{
  my ($self, $prefix) = @_;
  
  my $res = {tempfiles => $prefix} if $self->save_tempfiles;

  my $runfile = "$prefix.run";
  open FILE, $runfile or $self->throw("Can't open $runfile for reading");
  my $params = {};
  while (<FILE>){
    if (/^(.+)\s=\s(.+)$/){
      $params->{$1} = $2;
    }
  }
  close FILE;
  
  #keep run parameters.
  #$res->{params}=$params;
 
  #grab the temps from the dG file.
  open FILE, "$prefix.dG" or die "Can't open $prefix.dG for reading";
  my @temps;
  while (<FILE>){
    next if /^#T/;
    my ($temp) = /^(\d+).+$/;
    push @temps, $temp;
#    $res->{$temp}->{dG} = $rtlnz;
  }
  close FILE;

  #store these in the results, 
  #$res->{temps} = \@temps;


  my (%dH,%dG, %dS, %Tm, $img);
  foreach my $temp (@temps){
    
    #hybrid-ss-min generates ct files in the form prefix.temp.ct
    #unless there's only 1 temp, in which case its just prefix.ct
    #grrrr.
    my $ctfile = (scalar(@temps)> 1) ? "$prefix.$temp.ct" : "$prefix.ct" ;  
    
    #parse the dG,dH,dS,Tm from .ct as per melt.pl:
    #this needs to be DH for RNA and DHD for DNA.
#    my $suffix = 'DHD'; #energy rules. I have no idea, copied from melt.pl
    my $suffix = 'DH';
    $suffix = 'DHD' if (defined $params->{NA} and $params->{NA} eq 'DNA');

    open DG, "<$prefix.dG" or die $!;
    open DH, "ct-energy -s$suffix $ctfile |" or die $!;
    scalar <DG>;
    

    while (defined(my $dG = <DG>) and defined (my $dH = <DH>)) {
      $dG = (split /\s+/, $dG)[1];
      chomp $dH;
      my $dS = ($dH - $dG) / (273.15 + $temp);
      my $Tm = $dH / $dS;
      	
#      warn sprintf "%.1f\t%.1f\t%.1f\t%.1f\n", $dG, $dH, 1000 * $dS, $Tm - 273.15;
      
      #format and stick in res
      #$res->{dG}->{$temp} = sprintf ("%.1f", $dG);
      #$res->{dH}->{$temp} = sprintf ("%.1f", $dH);
      #$res->{dS}->{$temp} = sprintf ("%.1f", 1000 * $dS);
      #$res->{Tm}->{$temp} = sprintf ("%.1f", $Tm - 273.15);

      %dG = ($temp => sprintf ("%.1f", $dG));
      %dH = ($temp => sprintf ("%.1f", $dH));
      %dS = ($temp => sprintf ("%.1f", 1000 * $dS));
      %Tm = ($temp => sprintf ("%.1f", $Tm - 273.15));
    }		
		     
    close DH or die $!;
    close DG or die $!;

    
    #this requires Vienna tools - which should be converted 
    #into a perl wrapper and then tested for. Lets just leave
    #for now.
    #don't bother plotting unless we've got somewhere to put them
#    if ($self->plot_dir){
#      
#      #draw structure
#      system ("$convert_bin $ctfile > $prefix.$temp.vienna");
#      system ("$draw_bin   < $prefix.$temp.vienna");
#
#      #RNAplot writes to rna.ps and can't be convinced to write to stdout.
#      #use tempfile to get a name again
#      #probly naview is the way forward. This'll do for now.
#      my (undef,$plotfile) = $self->io->tempfile(DIR    => $self->plot_dir, 
#						 UNLINK => 0,
#						 SUFFIX => ".$temp.ps");
#
#      system ("mv rna.ps $plotfile") ;
#      $img = {$temp => "$plotfile"};
#      
#    }
  }


  $res =  Bio::Annotation::RNAFold::HybridSSMin->new(
	args   => $params,
	temps  => \@temps,
	dH     => \%dH,
	dG     => \%dG,
	dS     => \%dS,
	Tm     => \%Tm,
	prefix => $prefix,
       );
  
  $res->img($img) if $img;

  return $res;


}




sub _check_setup{
  my $self = shift;

  #does the exe exist?
  my $exe = $self->executable;
  $self->throw("$exe not found. Are program_dir and program_name set correctly?") unless -e $exe;


}

=head2 arguments

  Note       : Overriden from WrapperBase.
  Function   : Mostly a getter for constucting the hybrid-ss-min 
               command line argument string.
               Can also be used to set multiple command line arguments
               in one go
  Returns    : A string of command line arguments, like "--tmin=10 --maxbp=30" 
  Usage      : (set) $obj->arguments(tmin=>10, tmax=>40);
               (get) my $args = $obj->arguments;  

=cut

sub arguments{
  my ($self,%args) = @_;
  
  #try and store values if we've been given any.
  if (%args){
    foreach (keys %args){
      $self->$_($args{$_});
    }
  }
  
  #return concatenated values.
  my $args ='';
  foreach (@cmd_args){
    if ($self->$_){
      $args.="--$_";
      $args.='='.$self->$_.' ';
    }
  }
  return $args;
}


=head2 clear_arguments

   Function  : undef all the hybrid-ss-min arguments.
               Note that this only clears the hybrid-ss-min
               arguments, not everything in the object. Things
               like program_dir and seq_obj are still defined
   Usage     : $obj->clear_arguments;
   Returns   : true if successful.

=cut

sub clear_arguments{
  my $self = shift;
  
  foreach (@cmd_args){
     $self->$_(undef);
  }
  return 1;
}


#sub cleanup{
#  my $self = shift;
#  unless ($self->save_tempfiles){
#    if ($self->{_tmpdir}){
#      warn "Odd";
#      File::Path->rmtree( $self->{'_tmpdir'} ) or die "WTF:$!";
#    }
#  }
#}

1;
