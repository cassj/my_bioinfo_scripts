#!/usr/bin/perl
#biogeek bot for ##biogeeks, at least until we get a permanent channel.
#all it does (for now) it join the channel and keep it open
#code mostly thieved wholesale from http://poe.perl.org/?POE_Cookbook/IRC_Bots


use warnings;
use strict;
use POE;
use POE::Component::IRC;
sub CHANNEL () { "##biogeeks" }

# Create the component that will represent an IRC network.
my ($irc) = POE::Component::IRC->spawn();

# Create the bot session.  The new() call specifies the events the bot
# knows about and the functions that will handle those events.
POE::Session->create(
  inline_states => {
    _start     => \&bot_start,
    irc_001    => \&on_connect,
    irc_public => \&on_public,
  },
);

# The bot session has started.  Register this bot with the "magnet"
# IRC component.  Select a nickname.  Connect to a server.
sub bot_start {
  $irc->yield(register => "all");
  my $nick = 'biogeek_bot';
  $irc->yield(
    connect => {
      Nick     => $nick,
      Username => $nick,
      Ircname  => 'POE::Component::IRC biogeek bot',
      Server   => 'irc.freenode.net',
      Port     => '6667',
    }
  );
}

# The bot has successfully connected to a server.  Join a channel.
sub on_connect {
  $irc->yield(join => CHANNEL);
}

# The bot has received a public message.  Parse it for commands, and
# respond to interesting things.
#sub on_public {
#  my ($kernel, $who, $where, $msg) = @_[KERNEL, ARG0, ARG1, ARG2];
#  my $nick    = (split /!/, $who)[0];
#  my $channel = $where->[0];
#  my $ts      = scalar localtime;
#  print " [$ts] <$nick:$channel> $msg\n";
#  if (my ($rot13) = $msg =~ /^rot13 (.+)/) {
#    $rot13 =~ tr[a-zA-Z][n-za-mN-ZA-M];
#
#    # Send a response back to the server.
#    $irc->yield(privmsg => CHANNEL, $rot13);
#  }
#}

# Run the bot until it is done.

$poe_kernel->run();
exit 0;
