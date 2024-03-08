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
        my (%h,%AF);
        my $result = GetOptions(
                "file=s"         => \$opt{file},
		"help"  	 => \$opt{help}
        );
	if(!$result)            { die "ERROR: $! "}
        if($opt{help})          { print $HELP; exit 0 }

        ############################################################

        while(<>)
        {
		my @F=split /\t/;
		#next if(!/^#/ and (309<=$F[1] and $F[1]<=315 or 3100<=$F[1] and $F[1]<=3108));
		print;
        }
}
