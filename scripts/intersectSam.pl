#!/usr/bin/env perl 

use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that prints all the lines in a SAM file with queries present in a second file;

        EXAMPLE:
                $0 I.sam J.tab > O.sam

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
		"help"   => \$opt{help}
	);

	if(!$result)            { die "ERROR: $! "}
        if($opt{help})          { print $HELP; exit 0 }
        if(@ARGV<2)             { die "ERROR: insufficient number of arguments"}

        #########################################

        open(IN,$ARGV[1]) or die("ERROR: Cannot open input file".$!) ;
        while(<IN>)
        {
                next if(/^$/ or /^#/);

		chomp;
                my @F=split /\t/;
                $h{"$F[0]"}=1;
        }
	close(IN);
        last unless(%h);

        #########################################

        open(IN,$ARGV[0]) or die("ERROR: Cannot open input file".$!) ;

        while(<IN>)
        {
                if(/^@/) { print; next}

		chomp;
                my @F=split /\t/;
		print "$_\n" if($h{$F[0]});
        }
	close(IN);

	exit 0;
}


