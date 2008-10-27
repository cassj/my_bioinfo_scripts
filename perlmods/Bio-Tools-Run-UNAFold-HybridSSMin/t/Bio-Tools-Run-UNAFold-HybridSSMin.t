# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl Bio-Tools-Run-UNAFold-HybridSSMin.t'

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More qw/no_plan/;
use Test::Exception;
use Bio::Seq;
BEGIN { use_ok('Bio::Tools::Run::UNAFold::HybridSSMin') };

#########################

# Insert your test code below, the Test::More module is use()ed here so read
# its man page ( perldoc Test::More ) for help writing this test script.

# can we create an obj?
my $hsmin = new Bio::Tools::Run::UNAFold::HybridSSMin;
isa_ok ($hsmin, 'Bio::Tools::Run::UNAFold::HybridSSMin');

#check auto-created methods and validation
can_ok($hsmin, 'NA');
can_ok($hsmin,'tmin');
can_ok($hsmin,'tinc');
can_ok($hsmin, 'tmax');
can_ok($hsmin,'sodium');
can_ok($hsmin,'magnesium');
can_ok($hsmin,'polymer');
can_ok($hsmin,'suffix');
can_ok($hsmin,'output');
can_ok($hsmin,'prohibit');
can_ok($hsmin,'force');
can_ok($hsmin, 'energyOnly');
can_ok($hsmin, 'noisolate');
can_ok($hsmin, 'mfold');
can_ok($hsmin, 'maxbp');
can_ok($hsmin, 'constraints');
can_ok($hsmin, 'basepairs');
can_ok($hsmin, 'circular');
can_ok($hsmin, 'allpairs');
can_ok($hsmin, 'maxloop');
can_ok($hsmin, 'nodangle');
can_ok($hsmin, 'simple');
can_ok($hsmin, 'prefilter');
can_ok($hsmin, 'program_name');
can_ok($hsmin, 'program_dir');
can_ok($hsmin, 'seq_obj');


ok(!$hsmin->has_seq, 'has_seq returns false OK');

#Check defined sub-type checking:

#test PosInt checking 
lives_and (sub {$hsmin->maxbp(30), 30}, 'PosInt checking passes OK');
throws_ok {$hsmin->maxbp(0)} qr/ does not pass /, 'PosInt checking fails OK';

#test PosOrZeroInt checking
lives_and(sub {is $hsmin->tmin(10), 10 }, 'PosOrZeroInt checking passes OK' );
throws_ok {$hsmin->tmin(-3)} qr/ does not pass /, 'PosOrZeroInt checking fails negative OK';
throws_ok {$hsmin->tmin("rhubarb")} qr/ does not pass /, 'PosOrZeroInt checking fails string OK';
throws_ok {$hsmin->tmin(10.7)} qr/ does not pass /, 'PosOrZeroInt checking fails float OK';

#test ijk checking
lives_and( sub {is $hsmin->force('3,10,30'),'3,10,30'}, 'ijk checking passes OK');
throws_ok {$hsmin->force('tuesday')} qr/ does not pass /, 'ijk checking fails string OK';

#test PWMax checking
lives_and( sub {is $hsmin->mfold('10,100,70'), '10,100,70'}, 'PWMax checking passes OK');
throws_ok {$hsmin->mfold('monkey')} qr/ does not pass /, 'PWMax checking fails string OK';


#test ReadableFile checking
lives_and( sub {is $hsmin->constraints('t/Bio-Tools-Run-UNAFold-HybridSSMin.t'),'t/Bio-Tools-Run-UNAFold-HybridSSMin.t'},'ReadableFile checking passes OK');
throws_ok {$hsmin->constraints('curtains')} qr/ does not pass /, 'ReadableFile checking fails string OK';

#test Prefilter checking
lives_and ( sub {is $hsmin->prefilter('1,3'), '1,3'}, 'Prefilter checking passes OK');
throws_ok {$hsmin->prefilter('heffalump')} qr /does not pass/, 'Prefilter checking fails string OK';
throws_ok {$hsmin->prefilter('1,3,kitten')} qr/does not pass/, 'Prefilter checking fails mixed OK';


#test Dir checking
lives_and ( sub {is $hsmin->program_dir('.'), '.'}, 'Dir checking passes OK');
throws_ok {$hsmin->program_dir('nuts')} qr/does not pass/, 'Dir checking fails string OK';
throws_ok {$hsmin->program_dir('t/Bio-Tools-Run-UNAFold-HybridSSMin.t')} qr/does not pass/, 'Dir checking fails file OK';

#test SeqObj obj check.
my $seq = Bio::Seq->new();
throws_ok { $hsmin->seq_obj($seq) } qr/does not pass/, 'Bio::SeqI checking fails no seq OK';
$seq->seq('AGGGTGTCCCTACTAGCATCGATCATATATTGCTCGACTCTGGATATTCATCTGATGCATGCA');
$seq->display_name('test_sequence');
lives_and ( sub {is $hsmin->seq_obj($seq), $seq}, 'Bio::SeqI checking passes OK');
throws_ok { $hsmin->seq_obj('blah') } qr/does not pass/, 'Bio::SeqI checking passes OK';
ok($hsmin->has_seq, 'has_seq returns true OK');




#try setting some args with arguments.
$hsmin->arguments(seq_obj=>$seq, constraints=>undef);
$hsmin->tmin(35);
$hsmin->tmax(40);

#try clearing everything
$hsmin->clear_arguments; 

$hsmin->verbose(1);

$hsmin->tmin(35);
$hsmin->tmax(40);
$hsmin->mfold(1);
$hsmin->save_tempfiles(1);
$hsmin->plot_dir('/tmp');

#run the program as melt.pl
my $res = $hsmin->run;
$hsmin->cleanup;

# #try clearing everything and running with defaults
# $hsmin->clear_arguments;
# $hsmin->save_tempfiles(1);
# $res = $hsmin->run;
