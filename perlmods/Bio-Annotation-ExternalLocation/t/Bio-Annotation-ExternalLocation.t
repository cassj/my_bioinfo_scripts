# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl Bio-Annotation-ExternalLocation.t'

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More tests => 36;
BEGIN { use_ok('Bio::Annotation::ExternalLocation') };

#########################

# Insert your test code below, the Test::More module is use()ed here so read
# its man page ( perldoc Test::More ) for help writing this test script.

my $loc = Bio::Annotation::ExternalLocation->new();
isa_ok($loc, 'Bio::Annotation::ExternalLocation');

#RangeI
can_ok($loc, 'start');
can_ok($loc, 'end');
can_ok($loc, 'strand');
can_ok($loc, 'length');

#AnnotationI
can_ok($loc, 'as_text');
can_ok($loc, 'hash_tree');
can_ok($loc, 'tagname');

#Other Methods
can_ok($loc, 'authority');
can_ok($loc, 'taxon');
can_ok($loc, 'url');
can_ok($loc, 'coord_sys_type');
can_ok($loc, 'coord_sys_version');
can_ok($loc, 'coord_sys_id');


#Check setting of stuff at object creation
$loc = Bio::Annotation::ExternalLocation->new(
		'-start'             => 1,
                '-end'               => 100,
                '-strand'            => '-1',
                '-authority'         => 'NCBI',
                '-coord_sys_type'    => 'chromosome',
                '-coord_sys_version' => 37,
                '-coord_sys_id'      => 'X',
);

is($loc->start,1,'-start ok');
is($loc->end,100,'-end ok');
is($loc->strand, '-1', '-strand ok');
is($loc->authority,'NCBI', '-authority ok');
is($loc->coord_sys_type, 'chromosome', '-coord_sys_type ok');
is($loc->coord_sys_version, 37, '-coord_sys_version ok');
is($loc->coord_sys_id, 'X', '-coord_sys_id ok');

#And again without the bioperl dashes.
$loc = Bio::Annotation::ExternalLocation->new(
		start             => 1,
                end               => 100,
                strand            => '-1',
                authority         => 'NCBI',
                coord_sys_type    => 'chromosome',
                coord_sys_version => 37,
                coord_sys_id      => 'X',
);

is($loc->start,1,'start ok');
is($loc->end,100,'end ok');
is($loc->strand, '-1', 'strand ok');
is($loc->authority,'NCBI', 'authority ok');
is($loc->coord_sys_type, 'chromosome', 'coord_sys_type ok');
is($loc->coord_sys_version, 37, 'coord_sys_version ok');
is($loc->coord_sys_id, 'X', 'coord_sys_id ok');

#check tagname
is($loc->tagname,'ExternalLocation', 'default tagname ok');
$loc->tagname('Foo');
is($loc->tagname, 'Foo', 'setting tagname ok');

#test taxon if we can get a db connection
SKIP: {
  my $taxon;
  eval { 
    require Bio::DB::Taxonomy;
    my $db = new Bio::DB::Taxonomy(-source => 'entrez');
    my $taxonid = $db->get_taxonid('Homo sapiens');
    $taxon = $db->get_taxon(-taxonid => $taxonid);
  };

  skip "No DB connection for taxon retrieval", 4 unless $taxon;
  $loc->taxon($taxon);
  is($loc->taxon->id, $taxon->id, 'taxon ok');
  isa_ok($loc->taxon, 'Bio::Taxon', 'taxon class ok');

  #check length calc 
  is($loc->length, 100, 'length calculated ok');

  #check as_text
  is($loc->as_text, 'Auth: NCBI; Species: Homo sapiens; Version: 37; Type: chromosome; ID: X; 1 to 100 (100 bases on strand -1).', 'as text ok');

  #get hashes to compare
  my %test = (
	      'length'            => 100,
	      'coord_sys_version' => 37,
	      'coord_sys_type'    => 'chromosome',
	      'strand'            => '-1',
	      'species'           => 'Homo sapiens',
	      'coord_sys_id'      => 'X',
	      'authority'         => 'NCBI',
	      'end'               => 100,
	      'start'             => 1
	     );

  my %hash_tree = %{$loc->hash_tree};

  #sort hashes by keys and flatten them to a string
  my $test = join '', map {$_ , $test{$_}} sort {$a cmp $b} keys %test;
  my $hash_tree =  join '', map {$_, $hash_tree{$_}} sort {$a cmp $b} keys %hash_tree;

  #and compare them
  is($test, $hash_tree, 'hash_tree ok');
}

