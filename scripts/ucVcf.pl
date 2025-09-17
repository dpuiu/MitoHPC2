#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that ...

        EXAMPLE:
                cat I.vcf -min min_count | $0 

~;

###############################################################################
#
# Main program
#
###############################################################################


MAIN:
{
        my %opt;
        my $result = GetOptions(
                "help"          => \$opt{help},
        );
        if(!$result)            { die "ERROR: $! "}
        if($opt{help})          { print $HELP; exit 0 }

        ##########################################################################

	my (%header,%line,%count);
	while(<>)
	{
		if(/^#/)
		{
			print unless($header{$_}); 
			next;
		}

		my @F=split /\t/;
		($F[3],$F[4])=(uc($F[3]),uc($F[4]));
		print join "\t",@F;
	}
}

