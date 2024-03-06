#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that prints all the lines in a file with keys present in a second file;
  by default the keys are the values in the 1st columns

        EXAMPLE:
                $0 I.vcf J.vcf [ -i i -j j ]

~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
        my (%opt,%h);
	my($i,$j)=(0,0);

        # validate input parameters
        my $result = GetOptions(
		"i=i"	 => \$i,
		"j=i"	 => \$j,
		"header" => \$opt{header},
                "help"   => \$opt{help}
        );

	if(!$result)            { die "ERROR: $! "}
        if($opt{help})          { print $HELP; exit 0 }
        if(@ARGV<2)             { die "ERROR: insufficient number of arguments"}

        #########################################

        open(IN,$ARGV[1]) or die("ERROR: Cannot open input file".$!) ;
        while(<IN>)
        {
                next if(/^$/ or /^\@/ or /^#/);

		chomp;
                my @F=split /\t/;
                $h{"$F[$j]"}=1;
        }
	close(IN);
        last unless(%h);

        #########################################

        open(IN,$ARGV[0]) or die("ERROR: Cannot open input file".$!) ;

        while(<IN>)
        {
                if(/^\@/ or /^#/ or $.==1 and $opt{header}) { print; next}
		
		chomp;
                my @F=split /\t/;
		print "$_\n" if($h{$F[$i]});
        }
	close(IN);

	exit 0;
}


