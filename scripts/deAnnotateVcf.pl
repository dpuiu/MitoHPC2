#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that removes all the annotation (INFo tags) from a VCF file except for INDEL,GT,DP,AF,SM

        EXAMPLE:
                $0 I.vcf > O.vcf

~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
        my %options;
        my $result = GetOptions(
                "help"  =>	\$options{help}
	);

        if(!$result)            { die "ERROR: $! "}
        if($options{help})	{ print $HELP; exit 0 }

        #########################################

        my @tags=("INDEL","GT","DP","AF","SM");

        while(<>)
        {
                if(/^##INFO=<ID=(\w+),/)
		{
			my $keep=0;
			next unless $1 ~~ @tags;
			print;
		}
		elsif(/^#/)   { print }
		else
		{
			chomp;
	                my @F=split /\t/;
			$F[7]=";$F[7];";
			my $F7="";
			foreach  (@tags)
			{
				$F7.="$1;" if($F[7]=~/;($_);/ or $F[7]=~/;($_=.+?);/);
			}
			$F7="." unless($F7);
			$F7=~s/;$//;
			$F[7]=$F7;

			print join "\t",@F;print "\n";
		}
        }

	exit 0;
}


