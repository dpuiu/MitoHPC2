#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that filters columns 1-8 from a VCF files

        EXAMPLE:
                cat I.vcf  | $0 

~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
	# define variables
	my %opt;
	my $result = GetOptions(
		"help"          => \$opt{help}
        );
        if(!$result)            { die "ERROR: $! "}
        if($opt{help})          { print $HELP; exit 0 }

	#######################################################

	while(<>)
	{
		if(/^##/)
		{
			next if(/^##FILTER/ or /^##FORMAT/);
			print;
		}
		else
		{
			my @F=split /\t/;
			print join "\t",@F[0..7]; print "\n";
		}
	}
	exit 0;
}

