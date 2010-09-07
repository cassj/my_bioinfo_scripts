#!/usr/bin/perl

use strict;
use warnings;

$_[0] = [qw/m n o o r r p p p q/];
$_[1] = [qw/n p p q r r r s/];

use Data::Dumper;
warn Dumper( keys%{{map{$_,1} grep(${{map{$_=>1}@{$_[0]}}}{$_},@{$_[1]})}} ) ;

print "\n\n";


#all this is doing is $_[1] %in% $_[0]
#so if you've got 3 'r's in and the first one's got 2, you get three.
#which is not what you were asked for
#grep (${ { map{$_=>1} @{$_[0]}}  }{$_}, @{$_[1]})






die Dumper  { map{$_=>++} @{$_[0]}};
  

#grep (     ,@{$_[1]}) );
