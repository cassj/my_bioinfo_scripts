#!/usr/bin/perl
#biogeek bot for ##biogeeks, at least until we get a permanent channel.
#all it does (for now) it join the channel and keep it open
#code mostly thieved wholesale from POE::Component::IRC perldoc


use warnings;
use strict;
use POE;
use POE::Component::IRC;

my $channel = '##biogeeks'; 
my $server = 'irc.freenode.net';
my $nick = 'biogeek_bot';
my $ircname = 'London Biogeeks Bot';


# Create the component that will represent an IRC network.
my ($irc) = POE::Component::IRC->spawn(
    nick => $nick,
    ircname => $ircname,
    server => $server
    ) or die "New life blooms in spring\nJust like the giant panda\nBot failed to spawn ";


POE::Session->create(
    package_states => [
	main => [ qw(_default _start irc_001 irc_public) ],
    ],
    heap => { irc => $irc },
    );

$poe_kernel->run();


sub _start {
    my $heap = $_[HEAP];
    
    # retrieve our component's object from the heap where we stashed it
    my $irc = $heap->{irc};
    
    $irc->yield( register => 'all' );
    $irc->yield( connect => { } );
    return;
}

sub irc_001 {
    my $sender = $_[SENDER];
    
    # Since this is an irc_* event, we can get the component's object by
    # accessing the heap of the sender. Then we register and connect to the
    # specified server.
    my $irc = $sender->get_heap();
    
    print "Connected to ", $irc->server_name(), "\n";
    
    # we join our channels
    $irc->yield( join => $channel);
    return;
}

sub irc_public {
    my ($sender, $who, $where, $what) = @_[SENDER, ARG0 .. ARG2];
    my $nick = ( split /!/, $who )[0];
     my $channel = $where->[0];

     if ( my ($rot13) = $what =~ /^rot13 (.+)/ ) {
	$rot13 =~ tr[a-zA-Z][n-za-mN-ZA-M];
	$irc->yield( privmsg => $channel => "$nick: $rot13" );
    }
    return;
}

# We registered for all events, this will produce some debug info.
sub _default {
    my ($event, $args) = @_[ARG0 .. $#_];
    my @output = ( "$event: " );

    for my $arg (@$args) {
	if ( ref $arg eq 'ARRAY' ) {
	    push( @output, '[' . join(', ', @$arg ) . ']' );
	}
	else {
	    push ( @output, "'$arg'" );
	}
    }
    print join ' ', @output, "\n";
    return 0;
}
