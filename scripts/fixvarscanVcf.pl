#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that  ...  varscan VCF output

        EXAMPLE:
                 cat I.vcf | $0 > O.vcf

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
                "file=s"         => \$opt{file},
		"help"  	 => \$opt{help}
        );
	if(!$result)            { die "ERROR: $! "}
        if($opt{help})          { print $HELP; exit 0 }

        ############################################################

        while(<>)
        {
		if(/^#/)
		{
			print;
			next;
		}

		chomp;
		my @F=split /\t/;
		next if($F[3]=~/N/);
		print "$_\n";
        }
}
