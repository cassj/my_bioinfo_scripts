#!/usr/bin/perl

# note to self - this needs testing. 


# example:
# perl motif_tester.pl -i ENSMUSG00000018476 -s mouse --five_up 1000 --five_down 1000 -v \
# -l 'Homo sapiens' --matrix_id MA0061.1 -d /space/motifs/ --threshold '75%'
#


use strict;
use warnings;

use Data::Dumper; #remove when complete.
use List::Util qw(min max sum)  ;
use Storable;
use Graph;
use IO::File;

use Bio::Perl;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::AlignIO;
use Bio::EnsEMBL::Registry;
use Bio::Tools::Run::RepeatMasker;
use Bio::Graphics;
use Bio::Tools::Run::Alignment::Muscle;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::Annotation::SimpleValue;
use Bio::Annotation::Collection;
use Bio::SeqFeature::Generic;
use TFBS::DB::FlatFileDir;
use Bio::SeqFeature::Gene::GeneStructure;
use Bio::SeqFeature::Gene::Transcript;


use Bio::Graphics;

use Getopt::Long;
use Pod::Usage;

my ($verbose, $nomask, $help, $species, $identifier);

my $five_up=0;
my $three_down=0;
#a positive value of five_down trumps any setting for three_down
my $five_down=0;
#a positive value for three_up trumps any value of five_up
my $three_up=0;

#Note that the version of the DB used is whatever version of the API you're using.
my $ensembl_host = 'ensembldb.ensembl.org';
my $ensembl_user = 'anonymous';
my ($ensembl_pass, $ensembl_port);

my $msa = "Clustalw";
my $msa_out = "clustalw";

my $classification;
my @limit_species = ();
my $available_species;
my $database_dir;
my @matrix_IDs   = ();
my $conservation = 70;
my $threshold    = 80;
my $window       = 50;
my $report_file;
my $ortho_pad  = 250;

my $img_size = 800;

#only return motifs if they are totally conserved in all species
my $only_conserved = 0;

GetOptions(
	   'help|h'             => \$help,
	   'verbose|v'          => \$verbose,
	   'nomask|n'           => \$nomask,
	   'species|s=s'        => \$species,
	   'identifier|i=s'     => \$identifier,
	   'five_up=i'          => \$five_up,
	   'five_down=i'        => \$five_down,
	   'three_down=i'       => \$three_down,
	   'three_up=i'         => \$three_up,
	   'ensembl_host=s'     => \$ensembl_host,
	   'ensembl_user=s'     => \$ensembl_user,
	   'ensembl_pass=s'     => \$ensembl_pass,
	   'ensembl_port=i'     => \$ensembl_port,
	   'msa|m=s'            => \$msa,
	   'msa_out|o=s'        => \$msa_out,
	   'classification|c=s' => \$classification,
	   'limit_species|l=s'  => \@limit_species,
           'available_species|a'=> \$available_species,
	   'database=s'        => \$database_dir,
	   'matrix_id=s'       => \@matrix_IDs,
	   'conservation=f'    => \$conservation,
	   'threshold_score=s' => \$threshold,
	   'window_size=i'     => \$window,
	   'image_size=i'      => \$img_size,
           'report_file|r=s'   => \$report_file,
	   'ortho_pad=i'       => \$ortho_pad,
	   'only_conserved'    => \$only_conserved
);
if($help)  {
    pod2usage(-exitstatus=>0, -verbose=>2);
}
elsif (!($species and $identifier) )  {
    pod2usage(1);
}
elsif (!$database_dir) {
    pod2usage(1);
}
if($classification){
	die "Sorry - the --classification (-c) option isn't working yet.";
}

my %limit_species = map {$_=>1} @limit_species;
$report_file = "report_$identifier.txt" unless $report_file;

my $report_fh = new IO::File;
$report_fh->open("> $report_file") or die "can't open file $report_file";

print $report_fh "Report for $identifier \n";


#define connection to the EnsEMBL DB
my $REGISTRY = 'Bio::EnsEMBL::Registry';
$REGISTRY->load_registry_from_db(
                                 '-host' => $ensembl_host,
                                 '-user' => $ensembl_user,
				 '-pass' => $ensembl_pass,
				 '-port' => $ensembl_port
                                );


print $report_fh "Report for gene $identifier\n\n";
print $report_fh "Parameters:\n";
print $report_fh "\tSpecies\t$species\n";
print $report_fh "\tfive_up\t$five_up\n";
print $report_fh "\tfive_down\t$five_down\n";
print $report_fh "\tthree_up\t$three_up\n";
print $report_fh "\tthree_down\t$three_down\n";
print $report_fh "\tmsa\t$msa\n";
print $report_fh "\tlimit species\t@limit_species\n";
print $report_fh "\tdatabase_dir\t$database_dir\n";
print $report_fh "\tmatrix_IDs\t@matrix_IDs\n";
print $report_fh "\tthreshold\t$threshold\n";
print $report_fh "\n\n";


#identify the orthologs of this gene.
my $member_adaptor = $REGISTRY->get_adaptor('Multi', 'compara', 'Member');
my $member = $member_adaptor->fetch_by_source_stable_id('ENSEMBLGENE',$identifier);

print "Got Member\n" if $verbose;

my $homology_adaptor = $REGISTRY->get_adaptor('Multi', 'compara', 'Homology');
my $homologies = $homology_adaptor->fetch_all_by_Member($member);

print "Got Homology\n" if $verbose;

# create a data structure of orthologs like:
# { ENSID => [$gene, $slice, $seq], 
#   ENSID => [$gene, $slice, $seq]...}
my $data = {$identifier => fetch_data($species, $identifier, $five_up, $five_down, $three_up, $three_down, $report_fh)};

#print the sequence out for reference
my $outfa = Bio::SeqIO->new(-file => ">$identifier.fa" , '-format' => 'fasta');
$outfa->write_seq($data->{$identifier}->[2]);

print $report_fh "\n\nOrthologs:\n\n";
foreach my $homology (@$homologies) {

  #we're only interested in the 1-to-1 orthologs 
  next unless $homology->description eq "ortholog_one2one";

  foreach my $member_attribute (@{$homology->get_all_Member_Attribute}) {
    my ($member, $attribute) = @{$member_attribute};

    #first member is always the query gene.
    next if $member->stable_id eq $identifier;

    #fetch classification info
    my $taxon = $member->taxon;

    if($available_species){
	print $taxon->binomial."\n";	
	next;
    }
    
    if( scalar(@limit_species) > 0 ){
     next unless $limit_species{$taxon->binomial}; 
     warn $taxon->binomial if $verbose;
    }
    
    my $pad_five_up = $five_up ? $five_up + $ortho_pad : 0;
    my $pad_five_down = $five_down ? $five_down + $ortho_pad : 0;
    my $pad_three_up = $three_up ? $three_up + $ortho_pad : 0;
    my $pad_three_down = $five_down ? $five_down + $ortho_pad : 0;

    $data->{$member->stable_id} = fetch_data($taxon->name, $member->stable_id, $pad_five_up, $pad_five_down, $pad_three_up, $pad_three_down, $report_fh);
  }

}


# ok, do the msa.
my @seqs = map {$data->{$_}->[2]} keys %$data;

# I *think* the default params look okay.
my $factory = "Bio::Tools::Run::Alignment::$msa"->new();
my $aln = $factory->align(\@seqs);


#print the aln out for reference
my $out = Bio::AlignIO->newFh('-file' => ">$identifier".'.aln', '-format' => 'clustalw');
print $out $aln;




#load the PWMs in
@matrix_IDs = split (",", join(',',@matrix_IDs));
my $db = TFBS::DB::FlatFileDir->connect($database_dir);
my $matrixset;

#use all if non specified
unless (scalar @matrix_IDs) {
  $matrixset = $db->get_MatrixSet(-matrixtype=>"PWM");
}
else {
  $matrixset = $db->get_MatrixSet(-IDs        => \@matrix_IDs,
				  -matrixtype => "PWM");
}



my $seq = $data->{$identifier}->[2];
my $seqstring = $seq->seq;



#### draw the results.

my $panel = Bio::Graphics::Panel->new(
				      -length => $seq->length,
				      -width  => $img_size,
				      -pad_left  => 50,
				      -pad_right => 50,
				      -pad_top   => 10,
				      -pad_bottom => 10,
				      -image_class=>'GD', 
					 );


my $full_length = Bio::SeqFeature::Generic->new(
                                                -start => 1,
                                                -end   => $seq->length,
                                               );


my @transcript_feats = $data->{$identifier}->[2]->get_SeqFeatures();


#bottom strand? Draw length track first and gene underneath
if($data->{$identifier}->[0]->strand != 1){
    $panel->add_track($full_length,
		      -glyph   => 'arrow',
		      -tick    => 2,
		      -fgcolor => 'black',
		      -east    => 1,
		      );
}

foreach (@transcript_feats){
  $panel->add_track($_,
		    -glyph   => 'transcript2',
		    -label => 1,
                   -glyph       => 'transcript2',
                    -bgcolor     => 'orange',
                    -fgcolor     => 'black',
                    -font2color  => 'red',
                    -bump        =>  +1,
                    -height      =>  12,
		   );
}

#top strand? Draw gene first then length track
if($data->{$identifier}->[0]->strand == 1){
    $panel->add_track($full_length,
		      -glyph   => 'arrow',
		      -tick    => 2,
		      -fgcolor => 'black',
		      -east    => 1,
		      );
}

#figure out the transcript start sites relative to the gene start site?

#Calculate the conservation for each col of the alignment 
my $nseqs = $aln->no_sequences;
my @aln_seqs = map { [ split '', $aln->get_seq_by_pos($_)->seq ]} 1..$nseqs;

my $min_cons = 0;
my $max_cons = $nseqs;


my @alphabet = ('A','T','G','C');
my $consfeature =  Bio::SeqFeature::Generic->new(-start        => 1,
                                                 -end          => $seq->length,
                                                 -display_name => 'conservation',
);

#which cols in the alignment correspond to positions in our original sequence? 
#ie just ignore gaps for now.
my @cols = map {$aln->column_from_residue_number($identifier, $_) } 1..$seq->length;
my $consft = Bio::SeqFeature::Generic->new(
					   -start        => 1,
					   -end          => $seq->length,
					   -display_name => 'Conservation',
                                           );

my $seqpos = 1;
my @cons;
foreach my $alnpos (@cols){

    my @col =  map {$aln_seqs[$_-1][$alnpos-1] } 1..$aln->no_sequences;

    #calculate letter frequencies.
    my @count;
    foreach my $char (@alphabet){ push @count, scalar(grep /$char/, @col)}
     
    #this should eventually be IC for the col, but for now, just take the
    #best % conservation. 
	
    my @p = map {$_ / $nseqs} @count;
    push @cons, max(@p);

    $consft->add_SeqFeature(Bio::SeqFeature::Generic->new(
							  -start => $seqpos,
							  -end   => $seqpos,
							  -score => max (@p),
                                               ));

     $seqpos++;
} 

my @constrack = $panel->add_track($consft,
				  -label => 1,
				  -glyph        => 'heat_map',
				  );

#for each matrix, search against this string.
my %seen;
my $it = $matrixset->Iterator;
while(my $pwm = $it->next){

  #check we don't have ID dups
  my $id = $pwm->ID;
  $seen{$id} ? die "Duplicate matrix IDs not allowed" : $seen{$id}++;

  print $report_fh "\n\nSearching for matrix $id (".$pwm->name.")\n\n";  
  my $class = $pwm->class;
  $class =~s/\t/ /g;   #tabs in class confuse pwmsearch
  $pwm->class($class);
  my $siteset = $pwm->search_seq(-seqstring=>$seqstring, 
				-threshold=>$threshold);

  my $fh = new IO::File;
  if ($fh->open("> $identifier.$id.hits.gff")) {
    print $fh  $siteset->GFF;
    $fh->close;
    }
  
  #get an iterator over the hits
  my $it = $siteset->Iterator;

  #a feature to store my TFBS features
  my $fw_name = my $rev_name =  $pwm->name || $id;
  $fw_name .=" forward";
  $rev_name .=" reverse";
  my $ft_fw = Bio::SeqFeature::Generic->new(-start        => 1, 
					    -end          => $seq->length, 
					    -display_name => $fw_name,
					    -label        => 1,
					);
  
  my $ft_rev = Bio::SeqFeature::Generic->new(-start        => 1, 
					    -end          => $seq->length, 
					    -display_name => $rev_name,
					    -label        => 1,
					);
 
  my @scores;
  #and the hits as subfeatures
  while( my $site = $it->next){
    
    my @site_cons =  @cons[($site->start)-1..($site->end)-1];
    @site_cons = reverse @site_cons if ($site->strand != $data->{$identifier}->[0]->strand);
    my $conservation = join ' ', @site_cons;
    
    if ($only_conserved){
      next unless (sum(@site_cons) == scalar(@site_cons));
    }

    
    my $site_ft = Bio::SeqFeature::Generic->new(
						-start          => $site->start,
						-end            => $site->end,
						-strand         => $site->strand,
						-score          => $site->rel_score,
						-display_name   => $site->rel_score,
						-primary        => $site->rel_score,
						-label          => 1,
					       );
      
    if ($site->strand == 1) {
      $ft_fw->add_SeqFeature($site_ft);
    }else{
      $ft_rev->add_SeqFeature($site_ft);
    }
    
    push @scores, $site->rel_score ;
    
    #start relative to sequence start
    my $site_start = $site->start;
    my $site_end = $site->end;
    
    #site relative to sequence start
    print $report_fh "\tHit: Relative Score=".$site->rel_score."  On sequence at position $site_start to $site_end in the ".$site->strand. " direction \n";
    
    #site in genome co-ords
    my $slice = $data->{$identifier}->[1];
    
    my $site_genome_start = $slice->strand == 1 ? $slice->start + ($site->start - 1) : $slice->start + ($slice->length - $site->end);
    my $site_genome_end = $slice->strand == 1 ? $slice->start + ($site->end - 1) : $slice->start + ($slice->length - $site->start);
    my $site_genome_strand = $site->strand * $slice->strand;
    print $report_fh "\tGenome Co-ordinates: chr".$slice->seq_region_name.':'.$site_genome_start.'-'.$site_genome_end.' on strand '.$slice->strand;
    
    print $report_fh "\tSequence: ".$site->seq->seq."\n";

    print $report_fh "\tConservation: ".$conservation."\n";
    print $report_fh "\n";
  }
  
  #NOTE TO SELF - the relative score isn't really that helpful - it doesn't factor in how likely we are to see this 
  #sequence anyway. We really want some kind of IC based measurement - i really need a better understanding of how PWM 
  #searches work. 

  #add fw and bw tracks for this motif
  if ($ft_fw->get_SeqFeatures){
      my $track_fw =  $panel->add_track(
					-glyph     => 'graded_segments',
					-label     => 1,
					-bgcolor   => 'blue',
#					-min_score => min(@scores),
#					-max_score => max(@scores)
					);
      $track_fw->add_feature($ft_fw);
  }
  
  if ($ft_rev->get_SeqFeatures){
      my $track_rev =  $panel->add_track(
					 -glyph     => 'graded_segments',
					 -label     => 1,
					 -bgcolor   => 'blue',
#					 -min_score => min(@scores),
#					 -max_score => max(@scores)
					 );
      $track_rev->add_feature($ft_rev);
  }
  
}







open FILE, ">$identifier.hits.png" or die "can't open image file";
print FILE $panel->png;
close FILE;


$report_fh->close;










#### Supporting Functions





sub fetch_data {
    my ($species, $identifier, $five_up, $five_down, $three_up, $three_down, $report_fh) = @_;
    
    # get an Ensembl gene adap for the species in question
    my $gene_ad = $REGISTRY->get_adaptor(
					 $species,
					 'core',
					 'Gene',
                                      );
    
    # get an Ensembl slice adap for the species in question
    my $slice_ad = $REGISTRY->get_adaptor(
					  $species,
					  'core',
					  'Slice',
					 );
    
    unless (defined $gene_ad){
      die "No Gene Adaptor for $species";
    }
    unless (defined $slice_ad){
      die "No Slice Adaptor for $species";
    }

    my $gene = $gene_ad->fetch_by_stable_id($identifier);
    my @transcripts = @{$gene->get_all_Transcripts};

    print $report_fh "Gene $identifier ($species) Chr".$gene->slice->seq_region_name.':'.$gene->start.'-'.$gene->end.'('.$gene->strand.")\n";
    foreach(@transcripts){
	print $report_fh "\t".$_->display_id.' Chr:'.$gene->slice->seq_region_name.':'.$_->start.'-'.$_->end."\n";
    }

    #calculate start and end on the *top* strand
    my ($start, $end);
    if($gene->strand == 1){
      $start = $three_up ? $gene->end - $three_up : $gene->start - $five_up;
      $end = $five_down ? $gene->start + $five_down : $gene->end + $three_down;
    }else{
      $start = $five_down ? $gene->end - $five_down : $gene->start - $three_down;
      $end = $three_up ? $gene->start + $three_up : $gene->end + $five_up;
    }


    print $report_fh "\tGenome slice Chr".$gene->slice->seq_region_name.":$start".'-'."$end\n\n";

    #fetch the slice on the top strand
    my $slice = $slice_ad->fetch_by_region($gene->slice->coord_system->name, $gene->slice->seq_region_name, $start, $end);

    #invert the slice if gene is on the -ve strand
    #note that this just changes the ->strand and revcomps the seq.
    #->start and ->end are still given on top strand as per usual
    if ($gene->strand == -1){
      $slice = $slice->invert;
      print "Inverting slice to negative strand\n" if $verbose;
    }


    #and get the sequece from the slice
    my $seq = $nomask ? $slice->seq() : $slice->get_repeatmasked_seq()->seq;

    # stick seq in a Bio::Seq 
    my $bioseq = Bio::Seq->new(
			    -seq => $seq,
			    -id  => $gene->display_id,
			   );

    
    foreach (@transcripts){
      
      my $trsc_location = &relative_to_slice_seq($_->start, $_->end, $slice);

      #don't bother if the transcript doesn't actually overlap our sequence
      next unless $trsc_location->{overlap};

      my $transcript_ft = Bio::SeqFeature::Generic->new(
							-start  => max($trsc_location->{start}, 1),
							-end    => min($trsc_location->{end}, $bioseq->length),
							-display_name => $_->stable_id,
						       );
 
      my $location = &relative_to_slice_seq($_->start, $_->end, $slice);
      
      #coding_start and _end are relative to the 5' of the transcript
      #my $coding_region_start  = $_->coding_region_start;
      #my $coding_region_end = $_->coding_region_end; 
 
      #although the exons are given in the order they appear in the transcript.
      my @exons = @{$_->get_all_Exons};
      foreach (@exons){
	
	#get exon position relative to sequence
	my $exon_location  =  &relative_to_slice_seq($_->start, $_->end, $slice);
	next unless $exon_location->{start_contained};
	
	$transcript_ft->add_SeqFeature(
				       Bio::SeqFeature::Generic->new(
								     -start        => $exon_location->{start},
								     -end          => min($exon_location->{end}, length($seq)),
								     -display_name => $_->stable_id,
								    )
				      );
      }



	$bioseq->add_SeqFeature($transcript_ft);
   }

    return [$gene, $slice, $bioseq];
}


#takes a region (start and end pos on the top strand)
#and a slice and returns the co-ordiates relative to the slice sequence, 
#(on whichever strand the slice sequence is).
#returns a hash-ref of 
#start 
#end 
#overlap (true if any of the region is in the slice)
#start_contained (true if the start of the region is in the slice) 
#end_contained (true if the end of the region is in the slice)
sub relative_to_slice_seq{
  my ($start, $end, $slice) = @_;

  my $new_start = $slice->strand == 1 ? $start - $slice->start : $slice->end - $end;
  $new_start++;

  my $new_end = $slice->strand == 1 ? $end -$slice->start : $slice->end - $start;
  $new_end++;

  my $start_contained = ($new_start > 1) && ($new_start < $slice->length) ? 1:0;
  my $end_contained = ($new_end > 1) && ($new_end < $slice->length) ? 1:0;

  my $overlap = ($start_contained
		 || $end_contained
		 || ( ($new_start < 1) && ($new_end > $slice->length) )
		)? 1:0;

  return ({start=>$new_start, end=>$new_end, overlap=>$overlap, start_contained=>$start_contained, end_contained=>$end_contained});
  
}


__END__
    

=head1 NAME

fetch_and_align_orthologs.pl - Fetch Ensembl Orthologs for a gene and align them

=head1 SYNOPSIS

./fetch_and_align_orthologs.pl --species <SPECIES> --identifier <ENSEMBL GENE ID>

=head1 OPTIONS

=over 8


=item B<-h or --help>

Print help and exit.

=item B<-a or --available_species>

List the species for which orthologs exist and exit

=item B<-i or -identifier>  <identifier>

REQUIRED: Ensembl Gene ID for the gene for which you want 
to retrieve and align orthologs.

=item B<-s or --species>  <species>

REQUIRED: Species corresponding to the identifier

=item B<-c or --classification>  <classification>

OPTIONAL: Species must fall into this classification to be included 
in the MSA. By default, all orthologs are used.

=item B<-l or --limit_species> '<binomial name>'

OPTIONAL: By default, all orthologs are used for the MSA, if you define species using
this option (using their binomial name, eg -l 'Homo sapiens' -l 'Mus musculus') then only
the specified species will be used.


=item B<-n   or --nomask>

OPTIONAL: By default, sequences are repeatmasked.
This switch suppresses masking.

=item B<-v   or  --verbose> 

OPTIONAL: Turn on verbose reporting


=item B<--ensembl_host> <ensembl_host>

OPTIONAL: defaults to the EBI ensembl server

=item B<--ensembl_port> <ensembl_port>

OPTIONAL: defaults to the EBI ensembl server

=item B<--ensembl_user> <ensembl_user>

OPTIONAL: defaults to anonymous

=item B<--ensembl_pass> <ensembl_pass>

OPTIONAL: undefined by default

=item B<--five_up> <integer value>

Amount (in bases) of extra sequence to retrieve upstream of the 
start of the gene. The same distance will be retrieved for orthologs. 
Defaults to 0.

=item B<--three_down> <integer value>

Amount (in bases) of extra sequence to retrieve downstream of
the end of the gene. The same distance will be retrieved for orthologs.
Defaults to 0.

=item B<--five_down> <integer value>

Amount (in bases) of sequence to retrieve downstream of the
start of the gene. This overrides any value of three_down. 
Useful if you want to retrieve promoter / 5'UTR regions.

=item B<--three_up> <integer value>

Amount (in bases) of sequence to retrieve upstream of the 
end of the gene. This overrides and value of five_up.
Useful if you want to retrieve 3'UTR regions

NOTE: Defining three_up and five_down makes no sense and will result 
in an error.

=item B<--ortho_pad> <integer value>

OPTIONAL: defaults to 250
Values of five_up, five_down, three_up, three_down are also used to determine 
the region to retrieve around the start of each orthologous gene. To maximise 
the possibility of a good alignment to the search gene, you can get extra sequence
from the ortholgs by setting --ortho_pad to the amount of extra bases you would like.


=item B<-m or --msa> <multiple sequence alignment algorithm>

OPTIONAL: Which Multiple Sequence Alignment tool to use.
Defaults to ClustalW. 
Can be any Bio::Tools::Run::Alignment::$msa

=item B<-o or --msa_out> <output format>

OPTIONAL: Output format for the MSA.
Defaults to clustalw.
Can be anything Bio::AlignIO knows about.


=item B<-d   or --database>  <directory name>

REQUIRED: Name of the FlatFileDir database directory to 
use for retrieving matrices. 
A sample database directory examples/SAMPLE_FlatFileDir 
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

Given an ensembl gene ID fetch that gene and it's ensembl-defined orthologs.
Produce a multiple alignment of the sequences.

=cut

