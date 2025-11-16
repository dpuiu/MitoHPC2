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
	$opt{t}=0;
        my $result = GetOptions(
                "file=s"         => \$opt{file},
		"help"  	 => \$opt{help},
		"t=i"		 => \$opt{t}
        );
	if(!$result )           { die "ERROR: $! "}
        if($opt{help})          { print $HELP; exit 0 }

	$opt{t}/=100 if($opt{t}>=1);
        ############################################################

        while(<>)
        {
                if(/^#/)
                {
                        print;
                }
                else
                {
			chomp;
                        my @F=split /\t/;
                        ($F[-2] eq "GT:DP:AD:AF") or die "Error";
			$F[-1]=~/^(.+):(.+):(\d+),(\d+):(.+)$/;
			my $AF=int(100*$4/($3+$4)+.5)/100;

			next if($AF<$opt{t});
			$AF=1 if(1-$AF<$opt{t});
                        print join "\t",(@F[0..8],"$1:$2:$3,$4:$AF\n");	
                }
        }
}
