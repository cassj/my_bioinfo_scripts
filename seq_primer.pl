#!/usr/bin/perl


=head1 seq_primer.pl

  This script attempts to automate the Buckley lab
  primer design protocol for a defined sequence

=cut


use strict;
use warnings;

use Data::Dumper; #remove when complete.
use List::Util qw(min max sum)  ;
use Storable;

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

my $report_file = "primer_report.txt";
my $img_file = "primer_img.png";
my $img_size = 800;

my $ignore_masking = 1;
my $run_blast = 0;

print "Enter the species:\n";
my $species = <>;
chomp $species;

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


# print settings
print REPORT "Primer design to given sequence \n\n";
print REPORT "Species: $species\n";
print REPORT "TM Threshold: $tm_threshold\n";

print REPORT "\nPrimer3 Paramters:\n";
print REPORT Dumper \%primer3_params;
print REPORT "\nUNAFold Parameters\n";
print REPORT Dumper \%mfe_params;
print REPORT "\nRepeatMasker Parameters\n";
print REPORT Dumper \@rm_params;
print REPORT "\n\n";


print "Please enter your sequence:\n";
my $target_seq = <>;
chomp $target_seq;


#turn the target seq into a Bio::Seq object
$target_seq = Bio::Seq->new(
			    -seq        => $target_seq,
			    -primary_id => 'target_seq'
);


#run repeatmasker
my $rpt_masker = Bio::Tools::Run::RepeatMasker->new(@rm_params);
my @masked_feats;
my $masked_seq;
eval {
  @masked_feats = $rpt_masker->run($target_seq);
  $masked_seq = $rpt_masker->masked_seq;
  
  foreach (@masked_feats){
    
    print REPORT "Masked region: ".$_->start." to ".$_->end. ' length '.($_->end - $_->start)."\n";
    print REPORT $_->primary_tag."\n";
    print REPORT $target_seq->subseq($_->start,$_->end)."\n\n";
  }

};
if ($@) {
  $masked_seq = $target_seq;
}




#run primer3 on it.
my $primer3 = Bio::Tools::Run::Primer3->new(-path => "/usr/bin/primer3_core",
					    -outfile => "temp.out",
					     );
  
$primer3->add_targets(%primer3_params);
$primer3->add_targets(-seq=>$masked_seq);
my $primer3_res = $primer3->run;
  
my $good_primers = 0;
my %primers;
my $best_primer;

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
      
      #make amplicon seqfeat
      my $amp = Bio::SeqFeature::Generic->new(
					      -start      => $lp->start,
					      -end        => $rp->end,
					      -primary_id => $lp->start.'_'.$rp->end,
					      -display_name => $lp->start.'_'.$rp->end,
					     );
      
      #and add the primers to it as subfeatures
      $amp->add_SeqFeature($lp);
      $amp->add_SeqFeature($rp);
      
      #replace the annotated_sequence of our PrimedSeq with seqfeature
      $this_res->{annotated_sequence} = $amp;
      
      #store primer
      $primers{$amp->start.'.'.$amp->end} = $this_res;
	
      #keep a count of primers 
      $good_primers++;

      #spit out some blurb
      print REPORT "Primer pair at ".$amp->start.' to '. $amp->end."\n";
      print REPORT "\t".$this_res->{left_primer}->seq->seq."\n";
      print REPORT "\tLeft Tm: ".$this_res->mfe_left_primer->Tm->{60}."\n";
      print REPORT "\t".$this_res->{right_primer}->seq->seq."\n";
      print REPORT "\tRight Tm: ".$this_res->mfe_right_primer->Tm->{60}."\n";
      print REPORT "\t".$this_res->amplicon->seq."\n";
      print REPORT "\tAmplicon Tm: ".$this_res->mfe_amplicon->Tm->{60}."\n\n";
    }
    else{
      #give up and move on to the next one.
      next;
    }
  }
}

unless($good_primers){
  die "No primers pass Tm thresholds\n";
}



### DRAWING


my $panel = Bio::Graphics::Panel->new(
				      -length => $masked_seq->length,
				      -width  => $img_size,
                                      -pad_left  => 50,
                                      -pad_right => 50,
				      -image_class=>'GD', 
				     );

#target sequence
my $full_length = Bio::SeqFeature::Generic->new(
                                                -start => 1,
                                                -end   => $masked_seq->length,
                                               );


$panel->add_track($full_length,
                  -glyph   => 'arrow',
                  -tick    => 2,
                  -fgcolor => 'black',
                  -double  => 1,
                 );


#this draws all the primers underneath
my @primed_seqs = map {$_->annotated_sequence} values %primers;

my $track3 = $panel->add_track(\@primed_seqs,
			       -glyph => 'segments',
			       -bgcolor => 'blue',
			       -label => 1
			      );




open FILE, ">$img_file" or die "can't open image file";
print FILE $panel->png;
close FILE;


