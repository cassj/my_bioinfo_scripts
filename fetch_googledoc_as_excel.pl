#!/usr/bin/perl -

use strict;
use warnings;

use LWP::UserAgent;
use XML::Simple;
use Time::Local;
 
 

# Create browser and XML objects, and send a request for authentication
my $objUA = LWP::UserAgent->new;
my $objResponse = Fetch($objUA, "");


my $sheetid = "0AttSWcBIdgTCdC0wRjVYUWlmQ21YdG5NMm1HN1gxNUE";
#GET  http://spreadsheets.google.com/feeds/download/spreadsheets/Export?key=resource_id&exportFormat=format

my $objResponse = $objUA->post(
	$hConfig{URL},
	{
	accountType	=> $hConfig{AccountType},
	Email		=> $hConfig{UserName},
	Passwd		=> $hConfig{Password},
	service		=> $hConfig{Serivce},
	source		=> $hConfig{Source},
	"GData-Version" => $hConfig{APIVersion},
	}
); 
 
# Fail if the HTTP request didn't work
die "\nError: ", $objResponse->status_line unless $objResponse->is_success;
 


$objResponse = Fetch($objUA, "http://spreadsheets.google.com/feeds/list/MY_SHEET_KEY/od6/private/full");
my $objWorksheet = $objXML->XMLin($objResponse, ForceArray => 1);
 
# Open the file handle to create the output file
open (fhWRITE, ">$hConfig{OutputFile}") || die "Could not write to output file: $!\n";
 
# Put a header on the file
print fhWRITE GetFileChunk($hConfig{Header});
 
# For each row in the Google Docs sheet, print the date, hours, participants and notes
foreach my $sRow (@{$objWorksheet->{entry}}) {
 
	# Print out the rows where the date falls within the current week
        # For this to work you have to use the D/M/Y format in your date field
        # Depending what you name your spreadsheet columns, the gsx:date, etc
        # elements in the XML will change. You can use Data::Dumper to print the
        # XML to see what you're getting back if needed for debugging
	if (IsDateInPeriod($sRow->{'gsx:date'}[0])) {
		print fhWRITE "<tr>\n";
 
		print fhWRITE "<td>" . $sRow->{'gsx:date'}[0] . "</td>\n";
		print fhWRITE "<td class=\"hours\">" . $sRow->{'gsx:timehours'}[0] . "</td>\n";
		print fhWRITE "<td>" . $sRow->{'gsx:participants'}[0] . "</td>\n";
		print fhWRITE "<td>" . $sRow->{'gsx:notes'}[0] . "</td>\n";	
 
		print fhWRITE "</tr>\n";
 
		# Accumulate the total hours
		$iTotalHours += $sRow->{'gsx:timehours'}[0];
 
	} 
}
 
 
# Add the totals row
print fhWRITE "<tr class=\"totals\">\n";
print fhWRITE "<td colspan=\"2\" class=\"hours\">$iTotalHours</td>\n";
print fhWRITE "<td colspan=\"2\">Total Hours for Period</td>\n";
print fhWRITE "</tr>\n";
 
# Put a footer on the file
print fhWRITE GetFileChunk($hConfig{Footer});
 
close fhWRITE || warn "Could not write to output file: $!\n";
 
#------------------------------------------------------------------------------
# Subroutines
#------------------------------------------------------------------------------
 
# Extract the authorization token from Google's return string
sub ExtractAuth {
	# Split the input into lines, loop over and return the value for the 
	# one starting Auth=
   	for (split /\n/, shift) { 
   		return $1 if $_ =~ /^Auth=(.*)$/; 
   	}
   	return '';
 }
 
# Fetch a URL
sub Fetch {
	# Create the local variables and pull in the UA and URL
	my ($objUA, $sURL) = @_;
 
	# Grab the URL, but fail if you can't get the content
	my $objResponse = $objUA->get($sURL);
	die "Failed to fetch $sURL " . $objResponse->status_line if !$objResponse->is_success;
 
	# Return the result
	return $objResponse->content;
}
 
 
# Bring in an external file chunk and print it out 
sub GetFileChunk {
	# Pull the file name to print into a local variable
	my $sFile = shift;
	my $sFileChunk;
 
	# Whip through the file and fetch it into an array
	open(fhREAD, $sFile) || die "Could not open $sFile: $!\n";
 
	while (<fhREAD>) {
		# Local variables
		my $sLine = $_;
		my $sReplacementString;	
 
		# Fetch the date range
		if (/BILLING__PERIOD/) {
			my ($sStartDate, $sEndDate) = FetchDateRange();
			$sLine =~ s/BILLING__PERIOD/$sStartDate - $sEndDate/;
		}
 
		# Append this line
		$sFileChunk .= $sLine;
	}
 
	# Close the file handle
	close(fhREAD) || warn "Could not close $sFile: $!\n";
 
	# Return the result
	return($sFileChunk);
}
 
