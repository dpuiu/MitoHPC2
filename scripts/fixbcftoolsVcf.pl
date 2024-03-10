#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that adds the AF tag to the bcftools VCF output

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

                if(/^#/)
                {
                        print;
                        if(/^##INFO=<ID=DP,/)
                        {
                                print "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n";
                                print "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">\n";
                        }
                }
                else
                {
			chomp;
                        my @F=split /\t/;
                        next if($F[3]=~/N/);

                        my $DP=$1 if(/DP=(\d+)/);
                        my $AF=int(1000*($3+$4)/($1+$2+$3+$4)+.5)/1000 if(/DP4=(\d+),(\d+),(\d+),(\d+)/);
                        my $GT=$1 if(/.+\t(.+?):/);
                        print join "\t",(@F[0..7],"GT:DP:AF","$GT:$DP:$AF\n");
                }
        }
}
