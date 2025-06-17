#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that joins 2 TSV files by 1st column

        EXAMPLE:
                $0 I.tsv J.tsv > O.tsv

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
		"help"    =>	  \$opt{help}
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
		next unless(@F>=2);

                $h{"$F[0]"}=join "\t",@F[1..@F-1];
        }
	close(IN);

        #########################################

        open(IN,$ARGV[0]) or die("ERROR: Cannot open input file".$!) ;
        while(<IN>)
        {
                if(/^$/ or /^#/) { print; next}

		chomp; 
                my @F=split /\t/;
		$h{"$F[0]"}="." unless($h{"$F[0]"});
		print $_,"\t",$h{"$F[0]"},"\n";
        }

	exit 0;
}

