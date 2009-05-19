use Spreadsheet::ParseExcel;
use Net::Google::Calendar;
use DateTime;
use DateTime::Format::Natural;

print "Filename? ";
my $file = <>;
chomp $file;

my $parser   = Spreadsheet::ParseExcel->new();
my $workbook = $parser->Parse($file);

#we're only going to look at the first worksheet
my ($worksheet) = $workbook->worksheets();
if (!$worksheet){
  die "no worksheet found. bugger.";
}
my ( $row_min, $row_max ) = $worksheet->row_range();
my ( $col_min, $col_max ) = $worksheet->col_range();

my $cal = Net::Google::Calendar->new;

#it doesn't like googlemail.com usernames, but they map to
#gmail anyhow.
$cal->login('noel.buckley.lab@gmail.com', 'neur0sc1ence');

#delete all existing future lab meetings.
my @existing_events = $cal->get_events(q=>'Lab Meeting', futureevents=>'true');
foreach (@existing_events) {$cal->delete_entry($_);}

my $now = DateTime->now; 
my $year = $now->year;

#and add the new ones.
for my $row ( $row_min+1 .. $row_max ) {

  my ($day, $date, $time, $presenter, $candt) = map {$_ ? $_->value : ''}
    map {$worksheet->get_cell( $row, $_ ) || ''}
    ($col_min .. $col_max);

  next unless ($day && $time && $presenter);
  
  
  #this can parse the date, but not the time. fix at some point
  my $parser = DateTime::Format::Natural->new;
  my $start = $parser->parse_datetime("$date $year");
  my $end  = $parser->parse_datetime("$date $year");

  my ($t1_h, $t1_m, $t2_h, $t2_m) = $time =~ /^\s*(\d+)\s*:\s*(\d+)\s*-\s*(\d+)\s*:\s*(\d+)\s*$/;

  $start->set_hour($t1_h);
  $start->set_minute($t1_m);
  $end->set_hour($t2_h);
  $end->set_minute($t2_m);

  next unless $start >= $now;

  #enter into google calendar
  my $entry = Net::Google::Calendar::Entry->new();
  $entry->title("Lab Meeting: $presenter $candt");
  $entry->content("Lab Meeting: $presenter $candt");
  $entry->location('First Floor Meeting Room, James Black Centre, 125 Coldharbour Lane, London, England, SE5 9NU');
  $entry->transparency('transparent');
  $entry->status('confirmed');
  $entry->when($start, $end );

  $cal->add_entry($entry);
}

