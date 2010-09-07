#!/usr/bin/perl


=head1 design_tiled_primers.pl

  This script attempts to automate the Buckley lab
  tiling primer design protocol

  It could very probably do with some optimisation.

  Eventually, this will be converted to a proper
  catalyst app with a shiny interface, but at the 
  moment, we have nowhere to put such a beast,
  were it finished.

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
use Bio::Seq::PrimedSeq::Plus; #this is mine. and probly bollocks.
use Bio::Graphics;
use Bio::Tools::Run::Alignment::Blat;
use Bio::Tools::Run::RepeatMasker;


my $v = 1;

my $report_file = "tiled_primer_report.txt";
my $img_file = "tiled_primer_img.png";


print "Please enter species:\n";
my $species = <>;
chomp $species;

print "Please enter Ensembl Gene ID:\n";
my $identifier = <>;
chomp $identifier;

#my $identifier = 'ENSMUSG00000043969'; # mouse Emx2

#my $species = 'human';
#my $identifier = 'ENSG00000189056'; #human reelin
#my $identifier = 'ENSG00000077782'; #human fgfr1
#my $species = 'mouse';
#my $identifier='ENSMUSG00000037771 '; #mouse VIAAT

my $window_shift = 50;
my $window_size = 200;
my $min_amplicon = 90;
my $max_amplicon = 150;
my $min_inter_amplicon = 400;
my $max_inter_amplicon = 600;
my $upstream_of_tss= 1000;
my $downstream_of_tss = 1000;
my $tm_threshold = 63;


my $opt_inter_amplicon = 300;
#$min_inter_amplicon + ($max_inter_amplicon-$min_inter_amplicon)/2 ;



# set params as per Manu's protocol:
# Apparently these are already quite lenient.


my %primer3_params = (
		      PRIMER_OPT_GC_PERCENT     => 60,
		      PRIMER_MIN_GC             => 40,
		      PRIMER_MAX_GC             => 80,
		      PRIMER_PRODUCT_OPT_SIZE   => 120,
		      PRIMER_PRODUCT_SIZE_RANGE => '50-200',

		      PRIMER_OPT_SIZE           => 20,
		      PRIMER_MIN_SIZE           => 18,
		      PRIMER_MAX_MISPRIMING     => 12,
		      PRIMER_MIN_TM             => 57,
		      PRIMER_SELF_ANY           => 4,
		      PRIMER_GC_CLAMP           => 0,
		      PRIMER_NUM_NS_ACCEPTED    => 0,
		      PRIMER_OPT_TM             => 60,
		      PRIMER_MAX_POLY_X         => 5,
		      PRIMER_SALT_CONC          => 50,
#		      PRIMER_OUTSIDE_PENALTY    => 0,
		      PRIMER_MAX_TM             => 63,
		      PRIMER_SELF_END           => 3,
		      PRIMER_MAX_DIFF_TM        => 100,
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


#defaults are reasonable in PrimedSeqPlus
#you do need to specify an e value though
#
# blast really isn't what we want for this. we're
# looking for identical matches, surely?
# although actually, we only really care if the 3' e
# nd match is exact.
my %blast_params=(
  '-expect' => 1000,
);


my @rm_params=(-species  => $species,
	       -nolow    => 1,
	       -path    => "/usr/local/RepeatMasker",
	       -verbose => 1
	   );



##Open your report file and put some intro stuff in it

open REPORT, ">$report_file" 
  or die "Can't open file $report_file for writing";

print REPORT "Tiled Primer design for $identifier\n";

print REPORT "Settings:\n"
    ."Report file: $report_file\n"
    ."Image file: $img_file\n"
    ."Species: $species\n"
    ."Ensembl ID: $identifier\n"
    ."Min Amplicon Size: $min_amplicon"
    ."Max Amplicon Size: $max_amplicon"
    ."Min Inter-amplicon Distance: $min_inter_amplicon\n"
    ."Max Inter-amplicon Distance: $max_inter_amplicon\n"
    ."Optimal Inter-amplicon Distance: $opt_inter_amplicon\n"
    ."Upstream of TSS: $upstream_of_tss\n"
    ."Downstream of TSS: $downstream_of_tss\n"
    ."TM threshold: $tm_threshold\n";
    
print REPORT Dumper(\%primer3_params);
print REPORT Dumper(\%mfe_params);

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
print REPORT 'Gene: '.$gene->display_id.' '.$gene->description."retrieved from Ensembl $species database\n" ;
print REPORT "\tChromosome ".$gene->slice->seq_region_name;
print REPORT ' (start: '.$gene->start.' end: '.$gene->end;
print REPORT ' strand: '.$gene->strand.")\n";

my @exons = @{$gene->get_all_Exons};
print REPORT 'This gene has '.@exons." Exons\n";
my %exons;

for (my $i=0; $i<=$#exons; $i++){
  $exons{$exons[$i]->display_id} = $i;
  print REPORT $i."\t".$exons[$i]->display_id.': '
               .$exons[$i]->start.'-'
               .$exons[$i]->end."\n";
}

#get the transcripts and TSSs
my @transcripts = @{$gene->get_all_Transcripts};

print REPORT 'Gene has '.@transcripts. " transcripts\n";
my @tss;
foreach my $transcript (@transcripts){
  my $seq = $transcript->seq;
  print REPORT $transcript->display_id
               ."\t".$transcript->start
               ."-".$transcript->end
	       ."\n";
  my $tss = $gene->strand == 1 ? $transcript->start : $transcript->end;
  my $start_exon = $transcript->start_Exon;
  print REPORT "\t starts at chromosome position $tss in exon ".$exons{$start_exon->display_id}.' ('.$start_exon->display_id.')'."\n";
#  push @tss, $tss;
}


#just take the first TSS and a few kB around it.
my ($first_tss,$start,$end);

if ($gene->strand==1){
    $first_tss = $gene->start;
    $start = $first_tss - $upstream_of_tss;
    $end = $first_tss + $downstream_of_tss;
}else{
    $first_tss = $gene->end;
    #switch start and end around - needs to be +ve strand for slice
    $end = $first_tss + $upstream_of_tss;
    $start = $first_tss - $downstream_of_tss;
}

print REPORT "First TSS is at position $first_tss on chromosome ".$gene->slice->seq_region_name."\n";
print REPORT "Sequence retrieved from $start (TSS-3kb) to $end (TSS+5kb)\n";


my $slice = $slice_ad->fetch_by_region("chromosome", $gene->slice->seq_region_name, $start, $end);
my $target_sequence = $slice->seq;

# stick slice seq in a Bio::Seq
my $seq = Bio::Seq->new( -seq => $target_sequence,
 			 -id  => 'promoter_region_of:'.$gene->display_id,
		       );

#and reverse it because the slice is on the +ve strand
$seq = $seq->revcom unless $gene->strand == 1;


print REPORT 'sequence has length '.$seq->length;
print REPORT "\n\nTarget Sequence is:\n".$seq->seq."\n\n";


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

#### Get Primer3 results for a sliding window
# Check with mfold.
 
my $count = 0;         # for naming the primers
my $pos = 1;           # sequence is indexed from 1, not 0
$window_size--;        # window size is inclusive of the 1st base,   

print REPORT 
  "\n\nSliding Primer3 window of size $window_size, shift: $window_shift\n";


# Windows overlap, so we can get duplicates. 
my %seen;

WINDOW: while ($pos+$window_size < $masked_seq->length){
  
  $count++;
  warn "Processing window $count";

  #get the sequence corresponding to this window
  my $end_pos = $pos+$window_size;
  my $subseq  = $masked_seq->subseq($pos,$end_pos);

  print REPORT "\tWindow: $pos to $end_pos\n";
  print REPORT "\t$subseq\n";
  
  #make it into a Bio::Seq
  $subseq = Bio::Seq->new
    ( 
     -seq         => $subseq,
     -primary_id  => 'promoter_region_of:'.$gene->display_id."_$pos-$end_pos" ,
    );


  # get a new primer3 instance - it doesn't clean up from previous
  # runs properly, so if no primers are found, those from the previous
  # run are not overwritten and I get a 'can't place primer on seq' error.

  my $primer3 = Bio::Tools::Run::Primer3->new(-path => "/usr/bin/primer3_core",
					      -outfile => "temp.out",
					     );
  
  $primer3->add_targets(%primer3_params);
  
  #give it to the primer3 factory
  $primer3->add_targets(-seq=>$subseq);
  my $primer3_res = $primer3->run;
  
  my $good_primers = 0;
  #check we've got some results. 
  if ($primer3_res->number_of_results){
    
  PRIMER: while (my $this_res = $primer3_res->next_primer){
      
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
	  
	#this_res isa PrimedSeq. which isa SeqFeatureI
	my ($l_start, $l_len) = 
	  $this_res->{'left_primer'}->{'PRIMER_LEFT'} =~/(.+),(.+)/;
	my ($r_start,$r_len) = 
	  $this_res->{'right_primer'}->{'PRIMER_RIGHT'} =~ /(.+),(.+)/;

	#bioperl positions are indexed from 0, results of 
	#primer 3 are indexed from 1, so:
	$l_start++;
	$r_start++;
	
	
	my $st = $l_start+$pos-1;
	my $end = $r_start+$r_len+$pos-1;

	# set the position of the whole primed seq
	# relative to the main sequence
	# the amplicon isn't a seq feature, it's a seq
	# so it can't have location.
	$this_res->start($st);
	$this_res->end($end);
	$this_res->strand(1);
	my $this_num = $good_primers+1;
	$this_res->display_name("$count.$this_num");
	
	# set the positions of the primers relative to
	# main sequence.
	$this_res->{left_primer}->start($st);
	$this_res->{left_primer}->end($st+$l_len);
	$this_res->{right_primer}->start($r_start+$pos-1);
	$this_res->{right_primer}->end($end);


	#don't bother if we've seen it already.
	unless ( $seen{$st.'_'.$end} ){
	
	  #add the primed_seq to the seq as a seqfeature
	  $seq->add_SeqFeature($this_res);

	  $good_primers++;
	  print REPORT $this_res->display_name,"\n";
	  print REPORT "\t\tPrimer pair at ".$this_res->start.' to '. $this_res->end."\n";
	  print REPORT "\t\t\t".$this_res->{left_primer}->seq->seq."\n";
	  print REPORT "\t\t\tLeft Tm: ".$this_res->mfe_left_primer->Tm->{60}."\n";
	  print REPORT "\t\t\t".$this_res->{right_primer}->seq->seq."\n";
	  print REPORT "\t\t\tRight Tm: ".$this_res->mfe_right_primer->Tm->{60}."\n";
	  print REPORT "\t\t\t".$this_res->amplicon->seq."\n\n";
	  print REPORT "\t\t\tAmplicon Tm: ".$this_res->mfe_amplicon->Tm->{60}."\n";

	  #make a note of it
	  $seen{$st.'_'.$end}++;
	}
      }
      else{
	#give up and move on to the next one.
	next PRIMER;
      }

    }
    unless($good_primers){
      print REPORT "\t\tNo primers pass Tm thresholds in window\n";
    }
  }
  else{
    print REPORT "\t\tPrimer3 found no primers in window\n";
  }
  
  #shift window over and increase count.
  $pos=$pos+$window_shift;
}

### This is where we should come back to if we have ditched
### a primer pair







##############
#optimal tiling sequence
#seq should be a multi-primed seq, which could then have
#a method to do this

#these primers have all passed the Tm test.
my @primers = $seq->get_SeqFeatures;

#we don't want to bother blasting them unless we have to, 
#so construct a best-tile first, and if it fails a blast, 
#redo this step w/out offending pair and blast any extra

#this should be a blat, but I'm not sure how to do it online

#primer hash keyed by the name, <window_num>.<primer3_rank>
my %primers = map{$_->display_name => $_} @primers;

#seq is $seq
my @tm;
my (@i_pos, @j_pos, @start_amp, @end_amp);
my $temp = 60;  

#primer3 can run for mult temps, but we know it hasn't.
#duplicates were not stored, so we know positions are unique.

my @amps = sort {$a <=> $b} keys %primers;
my $alpha = scalar @amps;


foreach (@amps) {

  my $this_one = $primers{$_};

  #get the positions
  push @i_pos,$this_one->start;
  push @j_pos, $this_one->end;

  #get the tm values and use an ave of primers and amp
  push @tm, ($this_one->mfe_left_primer->Tm->{60}
  + $this_one->mfe_right_primer->Tm->{60}
  + $this_one->mfe_amplicon->Tm->{60})/3;

  #while we're at it, keep a note of any which are
  #in a valid start or end position and have an edge
  push @start_amp, $_ if $this_one->start < $max_inter_amplicon;
  push @end_amp, $_ if $this_one->end > ($seq->length - $max_inter_amplicon);

}


# use a graph to find the best set of tiles
my $g = Graph->new();

#each amplicon as a node, weighted by Tm.
for (my $a = 0; $a<=$#amps; $a++){
  $g->add_weighted_vertex("$amps[$a]", $tm[$a]);
}



#draw an edge between amplicons which are the 
#right distance apart

#thru end posns
for (my $ind_j = 0; $ind_j < $#j_pos; $ind_j++){
  
  #i & j not necessarily monotonic, so start
  #from 0 thru start posns until too far away

  my $min_i = $j_pos[$ind_j] + $min_inter_amplicon;
  my $max_i= $j_pos[$ind_j] + $max_inter_amplicon;

  # through the possible i positions, while the distance between i and
  # j is still shorter than the maximum:
   for (my $ind_i = 0;  
	$ind_i < $#i_pos && $i_pos[$ind_i] < $max_i; 
	$ind_i++){
     
          #for any distance bigger than the minimum
         if ($i_pos[$ind_i] > $min_i ){

	   #calculate a weight for this edge (0 at optimal, 
	   #rising to 1 at min and max inter_amplicon distance
	   my $d = $j_pos[$ind_j]-$i_pos[$ind_i]-1;
	   my $w = abs($d-$opt_inter_amplicon) / 100;

           #create an edge between the corresponding amplicon nodes
           $g->add_weighted_edge("$amps[$ind_j]", "$amps[$ind_i]", $w);


         }
   }
  
  
}



#not sure this will optimise the vertex weights
#might need to mod the edge weight to use the Tms.

#doesn't really matter though, all the tms are ok,

my @ok_paths;
foreach my $start_amp (@start_amp){
  foreach my $end_amp (@end_amp){
    my @sp = $g->SP_Dijkstra($start_amp,$end_amp);
    push @ok_paths, \@sp if @sp;
  }
}

die "No valid tiles found. Try increasing your inter-amplicon distance or your acceptable Tm." unless @ok_paths;

my @path;
my  $ave_edge_weight ;
foreach my $p (@ok_paths){

  #Each path is a list of vertices.
  my @p = @$p;

  #get the ave edge weight
  my @e_w;
  for (my $i = 0; $i<$#p; $i++){
    push @e_w, $g->get_edge_weight($p[$i], $p[$i+1])
  }
  my $this_ave_edge_weight =  (sum @e_w)/($#e_w+1); 
  
  #if its the lowest yet seen, keep it
  if  (!$ave_edge_weight || $this_ave_edge_weight < $ave_edge_weight){
    $ave_edge_weight = $this_ave_edge_weight;
    @path = @p;
  }

}


my @primed_seqs = $seq->get_SeqFeatures;
@primed_seqs = sort {$a->display_name > $b->display_name} @primed_seqs;

#get the primers corresponding to that path.
my @good_path = map {$primers{$_}} @path;


###
# draw the primers along the sequence?
###

my $panel = Bio::Graphics::Panel->new(
				      -length => $seq->length,
				      -width  => 1200,
                                      -pad_left  => 50,
                                      -pad_right => 50,
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



#this draws the masked regions in black
my $mask_track = $panel->add_track(\@masked_feats,
				   -glyph => 'generic',
				   bgcolor => 'black'
				  );



#and this draws the ve
my $track1 = $panel->add_track(\@good_path,
			       -glyph => 'segments',
			       -bgcolor => 'red',
			       -label => 1
			      );

#this draws the seq-feats of the seq feats.
my $track2 = $panel->add_track(\@primed_seqs,
			       -glyph => 'segments',
			       -bgcolor => 'blue',
			       -label => 1
			      );


open FILE, ">$img_file" or die "can't open image file";
print FILE $panel->png;
close FILE;


# now check they don't blast somewhere stupid. 
# if they do, ditch that primed_seq from the hash and do over.

#print their sequences out:
print REPORT "\n\n\n SELECTED AMPLICONS\n\n";








#I don't think we want to do this. BLAT would be better, they're 
#genomic so it should find them OK.
#warn "Blasting primers. This can take a very long time, please be patient";
foreach my $node (@good_path){

  #print out the amplicon sequence
  print REPORT "Amplicon Seq: ".$node->amplicon->seq."\n";
  print REPORT "Start: ".$node->start." End: ".$node->end."\n";
  print REPORT "Amplicon Mfold Tm: ".$node->mfe_amplicon->Tm->{60}."\n";
  print REPORT "Left primer Mfold Tm: ".$node->mfe_left_primer->Tm->{60}."\n";
  print REPORT "Right primer Mfold Tm: ".$node->mfe_right_primer->Tm->{60}."\n";

#  #blast the primers
#  my $left = $node->blast_left_primer(run=>1,%blast_params);
#  my $right = $node->blast_right_primer(run=>1,%blast_params);
#
#  print REPORT "Left primer BLAST results:\n\t".$node->left_blast_rid."\n";
#  print REPORT "Right primer BLAST results: \n\t".$node->right_blast_rid."\n";
#
#  warn "end of that loop";
# 
#  #check for problems (not implemented yet)
#  warn "CONFLICT WARNING: Both primers hit the same non-target site" if 
#    $node->blast_conflicts;
#
}

close REPORT or die "Can't close report file $report_file";
