#!/usr/bin/perl


=head1 genomic_primer_by_pos.pl

  Design primers to a specific region of the genome

=cut


use strict;
use warnings;

use Data::Dumper; #remove when complete.
use List::Util qw(min max sum)  ;
use Storable;
use Graph;

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
my $run_blast = 0;

print "Please enter species\n";
my $species = <>;
chomp $species;

print "Please enter the chromosome\n";
my $chr = <>;
chomp $chr;

print "Please enter start genome coordinate\n";
my $slice_start = <>;
chomp $slice_start;

print "Please enter end genome coordinate\n";
my $slice_end = <>;
chomp $slice_end;


print "Please enter amount of padding to use\n";
my $pad = <>;
chomp $pad;

my $report_file = "report_$chr:$slice_start-$slice_end.txt";
my $img_file = "img_$chr:$slice_start-$slice_end.png";
my $img_size = 800;


print "TM Threshold (default 63):\n";
my $tm_threshold = <>;
chomp $tm_threshold;
$tm_threshold = 63 unless $tm_threshold =~ /^\d+$/;

print "Min Amplicon Size (default 80):\n";
my $min_amplicon = <>;
chomp $min_amplicon;
$min_amplicon = 80 unless $min_amplicon =~ /^\d+$/;

print "Max Amplicon Size (default 170):\n";
my $max_amplicon = <>;
chomp $max_amplicon;
$max_amplicon = 170 unless $max_amplicon =~ /^\d+$/;


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
                      PRIMER_MAX_DIFF_TM        => 1,
                      PRIMER_MAX_SIZE           => 27,
                      PRIMER_NUM_RETURN         => 10
                     );



#Again, as per Manu's protocol
my %mfe_params = (
		  NA        => 'DNA',
		  tmin      => 60,
		  tmax      => 60,
		  sodium    => 0.05,
		  magnesium => 0.003,
);


my @rm_params=(-species  => $species,
	       -nolow    => 1,
	       -path    => "/usr/local/RepeatMasker",
	       -verbose => 1
	   );

##Open your report file and put some intro stuff in it

open REPORT, ">$report_file" 
  or die "Can't open file $report_file for writing";

print REPORT "Primer design for $chr:$slice_start-$slice_end\n";

### Get the sequence. ###

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db( 
 				 -host => 'ensembldb.ensembl.org',
 				 -user => 'anonymous',
 				);

# get an Ensembl slice adap for the species in question
my $slice_ad = $registry->get_adaptor(
				      $species,
				      'core',
				      'Slice',
				     );


# print some info about the gene
print REPORT "Species: $species";
print REPORT "\tChromosome $chr";
print REPORT " (start: $slice_start end: $slice_end pad: $pad)\n";


$slice_start = $slice_start-$pad;
$slice_end=$slice_end+$pad;
my $slice = $slice_ad->fetch_by_region("chromosome", $chr, $slice_start, $slice_end);
my $target_sequence = $slice->seq;
print REPORT "\n\nTarget Sequence is:\n".$target_sequence."\n\n";

# stick seq in a Bio::Seq
my $seq = Bio::Seq->new(
			-seq => $target_sequence,
			-id  => "$chr:$slice_start-$slice_end"
		       );

print REPORT 'sequence has length '.$seq->length;

# Get a masked_seq copy with repeat regions masked as Ns

#if this fails to find any masked sequences, it just dies
#so try-catch it.
my $rpt_masker = Bio::Tools::Run::RepeatMasker->new(@rm_params);
my @masked_feats;
my $masked_seq;
eval {
  @masked_feats = $rpt_masker->run($seq);
  $masked_seq = $rpt_masker->masked_seq;
  
  foreach (@masked_feats){
    
    print REPORT "Masked region: ".$_->start." to ".$_->end. ' length '.($_->end - $_->start)."\n";
    print REPORT $_->primary_tag."\n";
    print REPORT $seq->subseq($_->start,$_->end)."\n\n";
  }

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
  warn "Primer3 found ".$primer3_res->number_of_results.' primers';

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
      print REPORT "\t\tPrimer pair at ".$amp->start.' to '. $amp->end."\n";
      print REPORT "\t\t\t".$this_res->{left_primer}->seq->seq."\n";
      print REPORT "\t\t\tLeft Tm: ".$this_res->mfe_left_primer->Tm->{60}."\n";
      print REPORT "\t\t\t".$this_res->{right_primer}->seq->seq."\n";
      print REPORT "\t\t\tRight Tm: ".$this_res->mfe_right_primer->Tm->{60}."\n";
      print REPORT "\t\t\t".$this_res->amplicon->seq."\n\n";
      print REPORT "\t\t\tAmplicon Tm: ".$this_res->mfe_amplicon->Tm->{60}."\n";
      
    }#end if
  }#end while
  

  unless($good_primers){
    die "No primers pass Tm thresholds\n";
  }
  
} #end if
else{
  die "\t\tPrimer3 found no primers\n";
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

$panel->add_track($full_length,
                  -glyph   => 'arrow',
                  -tick    => 2,
                  -fgcolor => 'black',
                  -double  => 1,
                 );

#this draws the original sequence region
my $orig_seq = Bio::SeqFeature::Generic->new(
					     -start => 1+$pad,
					     -end   => $seq->length-$pad,
					     -display_name => 'Original Sequence Region'
					    );
$panel->add_track($orig_seq,
                  -glyph   => 'generic',
                  -bgcolor => 'green',
                  -height  => 20,
                  -label   => 1,
                  -description => 1);
 



#this draws the optimal primers
my $track1 = $panel->add_track($best_primer->annotated_sequence,
			       -glyph => 'segments',
			       -bgcolor => 'red',
			       -label => 1
			      );


#if we're not using masking info for primer design.
#at least see where it lies in relation to the primers.
#if ($ignore_masking){
if ($#masked_feats){
  #this draws the masked regions in black
  my $track2 = $panel->add_track(\@masked_feats,
				 -glyph => 'generic',
				 bgcolor => 'black'
				);
}

#this draws all the primers underneath
#why is it not drawing anything? I know I have primerss
my @primed_seqs = map {$_->annotated_sequence} values %primers;
my $track3 = $panel->add_track(\@primed_seqs,
			       -glyph => 'segments',
			       -bgcolor => 'blue',
			       -label => 1
			      );




open FILE, ">$img_file" or die "can't open image file";
print FILE $panel->png;
close FILE;


#print their sequences out:
print REPORT "\n\n\n SELECTED AMPLICONS\n\n";




#print out the amplicon sequence
print REPORT "Amplicon Seq: ".$best_primer->amplicon->seq."\n";
print REPORT "Start: ".$best_primer->annotated_sequence->start." End: ".$best_primer->annotated_sequence->end."\n";
print REPORT "Amplicon Mfold Tm: ".$best_primer->mfe_amplicon->Tm->{60}."\n";
print REPORT "Left primer sequence: ".$best_primer->{left_primer}->seq->seq."\n";
print REPORT "Right primer sequence: ".$best_primer->{right_primer}->seq->seq."\n";
print REPORT "Left primer Mfold Tm: ".$best_primer->mfe_left_primer->Tm->{60}."\n";
print REPORT "Left primer Mfold dG: ".$best_primer->mfe_left_primer->dG->{60}."\n";
print REPORT "Right primer Mfold Tm: ".$best_primer->mfe_right_primer->Tm->{60}."\n";
print REPORT "Right primer Mfold dG: ".$best_primer->mfe_right_primer->dG->{60}."\n";


close REPORT or die "Can't close report file $report_file";

