#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that colapses lines with the same key; 
1st column : the key
2nd column : values separated by ";"


        EXAMPLE:
                $0 I.tsv > O.tsv

~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
        my (%opt,%h);

        # validate input parameters
        my $result = GetOptions(
	);

	if(!$result)            { die "ERROR: $! "}
        if($opt{help})          { print $HELP; exit 0 }

        #########################################

	my %data;

	while (<>) {
	    chomp;
	    my ($id, $value) = split(/\t/, $_);
	    push @{ $data{$id} }, $value;
}

	foreach my $id (sort keys %data) {
	    my $collapsed = join(";", @{ $data{$id} });
            $collapsed = "." unless($collapsed);
	    print "$id\t$collapsed\n";
	}

	exit 0;
}

