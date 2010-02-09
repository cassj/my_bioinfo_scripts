#!/usr/bin/perl


=head1 cdna_primer.pl

  This script attempts to automate the Buckley lab
  primer design protocol

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

my $report_file = "cdna_primer_report.txt";
my $img_file = "cdna_primer_img.png";
my $img_size = 800;

my $ignore_masking = 1;
my $run_blast = 0;

print "Enter the species:\n";
my $species = <>;
chomp $species;


print "Enter the Ensembl Gene ID:\n";
my $identifier = <>;
chomp $identifier;


my $min_amplicon = 90;
my $max_amplicon = 200;
my $tm_threshold = 65;


# set params as per Manu's protocol:
# Apparently these are already quite lenient.

#my %primer3_params = (
#		      PRIMER_OPT_GC_PERCENT     => 60,
#		      PRIMER_MIN_GC             => 40,
#		      PRIMER_MAX_GC             => 80,
#		      PRIMER_PRODUCT_OPT_SIZE   => 120,
#		      PRIMER_PRODUCT_SIZE_RANGE => $min_amplicon.'-'.$max_amplicon,
#		      PRIMER_OPT_SIZE           => 20,
#		      PRIMER_MIN_SIZE           => 18,
#		      PRIMER_MAX_MISPRIMING     => 12,
#		      PRIMER_MIN_TM             => 57,
#		      PRIMER_SELF_ANY           => 4,
#		      PRIMER_GC_CLAMP           => 0,
#		      PRIMER_NUM_NS_ACCEPTED    => 0,
#		      PRIMER_OPT_TM             => 60,
#		      PRIMER_MAX_POLY_X         => 5,
#		      PRIMER_SALT_CONC          => 50,
#		      PRIMER_MAX_TM             => 63,
#		      PRIMER_SELF_END           => 3,
#		      PRIMER_MAX_DIFF_TM        => 100,
#		      PRIMER_MAX_SIZE           => 27,
#		      PRIMER_NUM_RETURN         => 5
#		     );
#



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
		      PRIMER_NUM_RETURN         => 20
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

print REPORT "cDNA Primer design for $identifier\n";



### Get the sequence. ###

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

# get an Ensembl slice adap for the species in question
my $slice_ad = $registry->get_adaptor(
				      $species,
				      'core',
				      'Slice',
				     );

# get the gene.
my $gene =$gene_ad->fetch_by_stable_id($identifier);

# print some info about the gene
print REPORT 'Gene: '.$gene->stable_id.' '.$gene->description."retrieved from Ensembl $species database\n" ;
print REPORT "\tChromosome ".$gene->slice->seq_region_name;
print REPORT ' (start: '.$gene->start.' end: '.$gene->end;
print REPORT ' strand: '.$gene->strand.")\n";

#and some to the screen
print "\n\n".$gene->stable_id.' on strand '.$gene->strand.' at position '.$gene->start.'-'.$gene->end."\n";

#get the transcripts and TSSs
my @transcripts = @{$gene->get_all_Transcripts};
print REPORT 'Gene has '.@transcripts. " transcripts\n";

#get exons for each transcript
my $n = scalar(@transcripts);

my %trsc_exons;
foreach my $trsc (@transcripts){
  $trsc_exons{$trsc->stable_id} = $trsc->get_all_Exons;
}


my $trsc;
if ($n>1){
  print "\nThere are multiple transcripts for this gene: \n";

  foreach my $t (@transcripts){
    print "\n".$t->stable_id;
    print " Known Coding" if $t->is_known;
    print "\nExons:\n ";
    foreach my $e ( @{$trsc_exons{$t->stable_id}} ){
      print "\t".$e->stable_id.' '.$e->start.'-'.$e->end."\n";
    }
  }
  print "\nPlease enter the ID of the transcript you would like to use: ";
  
  $trsc = <STDIN>;
  chomp($trsc);
  while (! $trsc_exons{$trsc}){
    print "$trsc doesn't seem to be a valid transcript ID. Please try again or stop the program with ctrl-c: ";
    $trsc = <STDIN>;
    chomp($trsc);
  }
}
else{
  my $t = $transcripts[0];
  $trsc = $t->stable_id;
  print "\n".$t->stable_id;
  print " Known Coding" if $t->is_known;
  print "\nExons:\n ";
  foreach my $e ( @{$trsc_exons{$t->stable_id}} ){
    print "\t".$e->stable_id.' '.$e->start.'-'.$e->end."\n";
  }
}

print "\nUsing $trsc \n";

my %exons;
foreach ( @{$trsc_exons{$trsc}} ){
  $exons{$_->stable_id} = $_;
#  print "\t".$_->stable_id.' ';
} 
#print "\n\n";

my $e_n = scalar (keys %exons);

my $target_seq;
my @target_exon_ids;
my $exon_bound = 1;

#single exon?
if ($e_n == 1){
  print "Only one exon, so we'll just use that sequence. You will be unable to detect genomic contamination with this primer pair\n";
  @target_exon_ids = scalar (keys %exons);
  $exon_bound = 0;
} 

##just 2?
#if ($e_n == 2){
#  print "Only two exons, so we'll aim for that exon-exon boundary";
#  #need to make sure these are the right way round.
#   @target_exon_ids = keys %exons;
#}

#more than 2:
if ($e_n >= 2){
  print "Which exon-exon boundary should I aim for?\n ";
  print "Please enter the ID of exon1: ";
  my $exon1 = <STDIN>;
  chomp($exon1);
  while (! $exons{$exon1} ){
    print "Invalid exon id. Try again, or stop the program with ctrl-c: ";
    $exon1 = <STDIN>;
    chomp($exon1); 
  }
  print "Please enter the ID of exon2: ";
  my $exon2 = <STDIN>;
  chomp($exon2);
  while (! $exons{$exon2} ){
    print "Invalid exon id. Try again, or stop the program with ctrl-c: ";
    $exon2 = <STDIN>;
    chomp($exon2); 
  }  
  @target_exon_ids = ($exon1,$exon2);
}

print "\n\nOK: using exon(s) ";
foreach (@target_exon_ids){print "$_ "}
print "\n\n";


#ok, concat the sequence of the exons into our target;
foreach (@target_exon_ids){
  $target_seq .= $exons{$_}->seq->seq;
}


#get the transcript object and its Bio::Seq
($trsc) = grep {$_->stable_id eq $trsc }@transcripts;
my $full_seq= $trsc->seq;

print REPORT "\nFull Transcript Sequence:\n".$full_seq->seq."\n";
print REPORT "\nTarget Sequence:\n".$target_seq."\n\n";

#and check that they haven't put an exon boundary that 
#isn't really in the transcript sequence:
die "That sequence doesn't exist in the transcript. Are you sure those exons are spliced next to each other?" unless $full_seq->seq =~ /$target_seq/;


#get the start and end positions of the target:
my $target_start = $exons{$target_exon_ids[0]}->start;
my $target_end   = $exons{$target_exon_ids[-1]}->end;


#turn the target seq into a Bio::Seq object
$target_seq = Bio::Seq->new(
			       -seq        => $target_seq,
			       -primary_id => "Primer target seq from ".$trsc->stable_id,
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
      print REPORT "\t".$this_res->amplicon->seq."\n\n";
      print REPORT "\tAmplicon Tm: ".$this_res->mfe_amplicon->Tm->{60}."\n";
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


#exon boundary
my $eb = $exons{$target_exon_ids[0]}->length;

$eb = Bio::SeqFeature::Generic->new(
				     -start => $eb,
				     -end => $eb
);

$panel->add_track($eb,
		  -glyph   => 'diamond',
		  -fgcolor => 'green',
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


