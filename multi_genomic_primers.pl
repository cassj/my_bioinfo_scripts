#!/usr/bin/perl


=head1 multi_genomic_primer.pl

  Design primers to a specific region of the genome, checking the 
  appropriate strand.

  Takes an input file which contains a list of chr and pos

=cut


use strict;
use warnings;

use Data::Dumper; 
use List::Util qw(min max sum)  ;
use Storable;
use Graph;
use IO::File;
use DateTime;


use Bio::Perl;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::AlignIO;
use Bio::EnsEMBL::Registry;
use Bio::Tools::Run::Primer3;
use Bio::Tools::Run::RepeatMasker;
use Bio::Seq::PrimedSeq::Plus;
use Bio::Graphics;

my $v = 1;


my $ignore_masking = 0;

my $input_file = shift @ARGV;
die "\nUsage:\n\tperl multi_genomic_primers.pl Inputfile\n\n" unless $input_file;
chomp($input_file);

print "Please enter species\n";
my $species = <>;
chomp $species;

my $min_amplicon = 100;
my $max_amplicon = 150;
my $window = 100;
my $tm_threshold =   63;


#Diogo's strict parameters
my %primer3_params = (
		      PRIMER_OPT_GC_PERCENT     => 60,
		      PRIMER_MIN_GC             => 30,
		      PRIMER_MAX_GC             => 80,
		      PRIMER_PRODUCT_OPT_SIZE   => 120,
		      PRIMER_PRODUCT_SIZE_RANGE => "$min_amplicon - $max_amplicon",
		      PRIMER_OPT_SIZE           => 20,
		      PRIMER_MIN_SIZE           => 18,
		      PRIMER_MAX_MISPRIMING     => 12,
		      PRIMER_MIN_TM             => 58,
		      PRIMER_SELF_ANY           => 4,
		      PRIMER_GC_CLAMP           => 0,
		      PRIMER_NUM_NS_ACCEPTED    => 0,
		      PRIMER_OPT_TM             => 59,
		      PRIMER_MAX_POLY_X         => 5,
		      PRIMER_SALT_CONC          => 50,
		      PRIMER_MAX_TM             => 60,
		      PRIMER_SELF_END           => 2,
		      PRIMER_MAX_DIFF_TM        => 2,
		      PRIMER_MAX_SIZE           => 27,
		      PRIMER_NUM_RETURN         => 5
		     );


#as per Manu's protocol
my %mfe_params = (
		  NA        => 'DNA',
		  tmin      => 60,
		  tmax      => 60,
		  sodium    => 0.05,
		  magnesium => 0.003,
);


my @rm_params =(-species  => $species,
		-nolow    => 1,
		-path    => "/usr/local/RepeatMasker",
		-verbose => 1
	       );


my $img_size = 800;




#store the settings 

my $param_fh = new IO::File;
$param_fh->open("> parameters.txt") or die "Can't open file parameters.txt for writing";

print $param_fh "***** Primer 3 *****\n\n",Dumper (\%primer3_params),"\n\n";
print $param_fh "***** Unafold *****\n\n",Dumper (\%mfe_params),"\n\n";
print $param_fh "***** Repeat Masker *****\n\n",Dumper (\@rm_params),"\n\n";
print $param_fh "***** Other *****\n\n";
print $param_fh "Input File: $input_file";

$param_fh->close;

### Database setup
my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db( 
 				 -host => 'ensembldb.ensembl.org',
 				 -user => 'anonymous',
 				);

# get an Ensembl gene adap for the species in question
my $gene_ad = $registry->get_adaptor(
				     $species,
				     'core', 
				     'Gene',
				    );


# get an Ensembl gene adap for the species in question
my $trsc_ad = $registry->get_adaptor(
				     $species,
				     'core', 
				     'Transcript',
				    );

# get an Ensembl slice adap for the species in question
my $slice_ad = $registry->get_adaptor(
				      $species,
				      'core',
				      'Slice',
				     );



my $input_fh = new IO::File;
$input_fh->open("< $input_file") or die "Can't open file $input_file for reading";

my $start_time = DateTime->now;


#skip the header
<$input_fh>;
while(my $this_row = <$input_fh>){

  #ditch any quotes.
  $this_row =~ s/\"//g;
  my ($identifier, $chr, $genome_start, $genome_end, $strand, $pval, $probe_start, $probe_end, $genome_mid, $feature_start, $feature_end) = split("\t",$this_row);
  $chr =~s/chr//;
  
  warn "Processing for $identifier";
  my $pos = $genome_mid;
  
  my $slice_start = $pos-$window;
  my $slice_end = $pos+$window;
  
  my $report_file="Primer_Report_$identifier.txt";
  my $img_file = "Primer_Image_$identifier.png";
  
  my $report_fh = new IO::File;
  $report_fh->open(">$report_file") or die "Can't open report file $report_file for writing";
  print $report_fh "Primer design for $identifier\n";
  
#  my $gene = $gene_ad->fetch_by_stable_id($identifier);
#  warn "$identifier has conflicting strand info" if $gene->strand != $strand; 
#  
#  # print some info about the gene
#  print $report_fh 'Gene: '.$gene->display_id.' '.$gene->description."retrieved from Ensembl $species database\n" ;
#  print $report_fh "\tChromosome ".$gene->slice->seq_region_name;
#  print $report_fh ' (start: '.$gene->start.' end: '.$gene->end;
#  print $report_fh ' strand: '.$gene->strand.")\n";
#  
#  my @exons = @{$gene->get_all_Exons};
#  print $report_fh "Exons:\n";
#  for (my $i=0; $i<=$#exons; $i++){
#    print $report_fh $exons[$i]->display_id.' : ';
#    print $report_fh $exons[$i]->start.' - '. $exons[$i]->end;
#    print $report_fh "\n";
#  }
#  
#  print $report_fh "\n\nTranscripts";
#  my @transcripts = @{$gene->get_all_Transcripts};
#  if ($#transcripts){
#    foreach my $ts (@transcripts){
#      print $report_fh $ts->display_id.' : ';
#      my @texons = @{$ts->get_all_Exons};
#      print $report_fh join ' ', map {$_->stable_id} @texons;
#      print $report_fh "\n";
#    }
#    print $report_fh "\n\n";
#  }
  
  #ok, our slice actually needs to be in a specifc region, so:
  
  my $slice = $slice_ad->fetch_by_region("chromosome", $chr, $slice_start, $slice_end);
  
  my $target_sequence = $slice->seq;
  
  print $report_fh "Searching region: ".$slice->start.'-'.$slice->end."\n";
  print $report_fh "Binding Peak is at".$pos ."\n";
  print $report_fh "Target sequence saved as last_target.fa\n";
  print $report_fh "\n\nTarget Sequence is:\n".$target_sequence."\n\n";

  # stick seq in a Bio::Seq
  my $seq = Bio::Seq->new(
			  -seq => $target_sequence,
			  -id  => "Primer_for_$identifier",
			 );
    
  print $report_fh 'sequence has length '.$seq->length;
  
  # Get a masked_seq copy with repeat regions masked as Ns
  #if this fails to find any masked sequences, it just dies
  #so try-catch it.
  my $rpt_masker = Bio::Tools::Run::RepeatMasker->new(@rm_params);
  my @masked_feats;
  my $masked_seq;
  eval {
    @masked_feats = $rpt_masker->run($seq);
    $masked_seq = $rpt_masker->masked_seq;
  };
  if ($@) {
    $masked_seq = $seq;
  }
  
  
  
  #### Get PrimedSeqs using Primer3 ###

  # get a primer3 instance and check the exe exists
  my $primer3 = Bio::Tools::Run::Primer3->new(-path => "/usr/bin/primer3_core",
					      -outfile => "temp.out",
					     );
  unless ($primer3->executable) {
      die "primer3 can not be found. Is it installed?";
    }

  $primer3->add_targets(%primer3_params);
  
  ## avoid duplicates
  my %primers;
  
  #give it to the primer3 factory
  $primer3->add_targets(-seq=>$masked_seq);
  my $primer3_res = $primer3->run;
  
  my $good_primers = 0;
  
  #check we've got some results. 
  if ($primer3_res->number_of_results){
    while (my $this_res = $primer3_res->next_primer){
      
      #make it into a PrimedSeqPlus
      $this_res = Bio::Seq::PrimedSeq::Plus->coerce($this_res);
      $this_res->species($species);
      
      #run mfe on the primers and amplicon
      #note that the right primer sequence is revcomp of the
      #target sequence ie. 5' to 3' on the opp strand.
      $this_res->mfe_left_primer(%mfe_params);
      $this_res->mfe_right_primer(%mfe_params);
      $this_res->mfe_amplicon(%mfe_params);
      
      #check all Tms pass a threshold
      if ($this_res->mfe_tms_under(threshold=>$tm_threshold)){
	
	my ($lp,$rp) = $this_res->annotated_sequence->get_all_SeqFeatures;
	  
	#we'll also create a seqfeature for the amplicon
	#just name it with the position
	my $amp = Bio::SeqFeature::Generic->new(
						-start      => $lp->start,
						-end        => $rp->end,
						  -primary_id => $lp->start.'_'.$rp->end,
						-display_name => $lp->start.'_'.$rp->end,
					       );
	
	#and add the primers to it as subfeatures
	$amp->add_SeqFeature($lp);
	  $amp->add_SeqFeature($rp);
	
	#replace the annotated_sequence of our PrimedSeq 
	#with out amplicon.
	#this isn't a setter method, so cheat.
	$this_res->{annotated_sequence} = $amp;
	
	#if not, add it to the list, keyed by amplicon position
	#note this will overwrite any duplicates 
	$primers{$amp->start.'.'.$amp->end} = $this_res;
	
	#keep a count of primers in this window
	$good_primers++;
	
	#spit out some blurb
	print $report_fh "\t\tPrimer pair at ".$amp->start.' to '. $amp->end."\n";
	print $report_fh "\t\t\t".$this_res->{left_primer}->seq->seq."\n";
	print $report_fh "\t\t\tLeft Tm: ".$this_res->mfe_left_primer->Tm->{60}."\n";
	print $report_fh "\t\t\t".$this_res->{right_primer}->seq->seq."\n";
	print $report_fh "\t\t\tRight Tm: ".$this_res->mfe_right_primer->Tm->{60}."\n";
	print $report_fh "\t\t\t".$this_res->amplicon->seq."\n\n";
	print $report_fh "\t\t\tAmplicon Tm: ".$this_res->mfe_amplicon->Tm->{60}."\n";
	
      }#end if
    }#end while
    
    
    unless($good_primers){
      $report_fh->close;
      unlink($report_file) or warn "Couldn't unlink $report_file";
      #die "No primers pass Tm thresholds\n";
      next;
      }
    
  } #end if
  else{
      $report_fh->close;
      unlink($report_file) or warn "Couldn't unlink $report_file";
      next;
      #die "\t\tPrimer3 found no primers\n";
    }#end else
    
  
  
  
  ####
  # Pick the best primer pair.
    
  my $best_primer;
  my $temp = 60;  
    
  #primer3 can run for mult temps, but we know it hasn't.
    #duplicates were not stored, so we know positions are unique.
  
  my @amps = sort {$a <=> $b} keys %primers;
  my $alpha = scalar @amps;
  
    foreach (@amps) {
      
      my $this_one = $primers{$_};
      
      #get the tm values and use total of primers and amp
      my $tm = ($this_one->mfe_left_primer->Tm->{60}
		+ $this_one->mfe_right_primer->Tm->{60}
		+ $this_one->mfe_amplicon->Tm->{60});
      
      if (! defined($best_primer) || 
	  $tm <
	  ($best_primer->mfe_left_primer->Tm->{60}
	   + $best_primer->mfe_right_primer->Tm->{60}
	   + $best_primer->mfe_amplicon->Tm->{60}) 
	 ){
	$best_primer = $this_one;
      }
    }




    ####
    ## draw the primers along the sequence?
    ###
    
    my $panel = Bio::Graphics::Panel->new(
					  -length     => $seq->length,
					  -width      => $img_size,
					  -pad_left   => 10,
					  -pad_right  => 10,
					  -pad_top    => 10,
					  -pad_bottom =>10,
					  -image_class=>'GD', 
					 );
    
    my $full_length = Bio::SeqFeature::Generic->new(
						    -start => 1,
						    -end   => $seq->length,
						   );
    
    
    #at this stage, pos is still relative to the genomic +1 seq.
    #shift to 1..$seq->length frame
    $pos = $strand == 1 ?
      $pos - $slice->start + 1:
      $slice->end - $pos + 1;
    

    
    $pos = Bio::SeqFeature::Generic->new(
					 -start => $pos,
					 -end => $pos
					);
    
    $panel->add_track($pos,
		      -glyph   => 'diamond',
		      -fgcolor => 'green',
		      -label   => 'BS Peak'
		     );

    
    
    $panel->add_track($full_length,
		      -glyph   => 'arrow',
		      -tick    => 2,
		      -fgcolor => 'black',
		      -double  => 1,
		     );

    
    #this draws the optimal primers
    my $track1 = $panel->add_track($best_primer->annotated_sequence,
				   -glyph => 'segments',
				   -bgcolor => 'red',
				   -label => 1
				  );


    if ($#masked_feats){
      #this draws the masked regions in black
      my $track2 = $panel->add_track(\@masked_feats,
				     -glyph => 'generic',
				     bgcolor => 'black'
				    );
    }

    #this draws all the primers underneath
    my @primed_seqs = map {$_->annotated_sequence} values %primers;
    my $track3 = $panel->add_track(\@primed_seqs,
				   -glyph => 'segments',
				   -bgcolor => 'blue',
				   -label => 1
				  );


    my $img_fh = new IO::File;
    $img_fh->open(">$img_file") or die "Can't open report file $img_file for writing";

    
    print $img_fh $panel->png;
    $img_fh->close;
    


    print $report_fh "\n\n\n SELECTED AMPLICON\n\n";
    print $report_fh "Amplicon Seq: ".$best_primer->amplicon->seq."\n";
    print $report_fh "Start: ".$best_primer->annotated_sequence->start." End: ".$best_primer->annotated_sequence->end."\n";
    print $report_fh "Amplicon Mfold Tm: ".$best_primer->mfe_amplicon->Tm->{60}."\n";
    print $report_fh "Left primer sequence: ".$best_primer->{left_primer}->seq->seq."\n";
    print $report_fh "Right primer sequence: ".$best_primer->{right_primer}->seq->seq."\n";
    print $report_fh "Left primer Mfold Tm: ".$best_primer->mfe_left_primer->Tm->{60}."\n";
    print $report_fh "Left primer Mfold dG: ".$best_primer->mfe_left_primer->dG->{60}."\n";
    print $report_fh "Right primer Mfold Tm: ".$best_primer->mfe_right_primer->Tm->{60}."\n";
    print $report_fh "Right primer Mfold dG: ".$best_primer->mfe_right_primer->dG->{60}."\n";

  $report_fh->close;


}

$input_fh->close;


my $end_time = DateTime->now;
my $duration = $end_time-$start_time;

warn "Processing time: ".$duration->hours.':'.$duration->minutes.':'.$duration->seconds;
