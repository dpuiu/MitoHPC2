#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that converts between VCF (con)catenated formats

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
        # define variables
        my %opt;

        my $result = GetOptions(
                "help"  =>	\$opt{help}
                );

        if(!$result)    { die "ERROR: $! "}
        if($opt{help})  { print $HELP; exit 0 }

	my %GT;
	while(<>)
	{
		if(/^#/)
		{
			s/INFO/FORMAT/ if(/^##INFO=<ID=DP,/ or /##INFO=<ID=AF,/);
			s/FORMAT/INFO/ if(/^##FORMAT=<ID=SM,/); 

			print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" if(/^##FORMAT=<ID=DP,/);
			print;
			next;
		}
		
		chomp;
		my @F=split /\t/;
		if($F[8] eq "SM")
		{
			my ($GT,$DP,$AD,$AF,$ANNOTATION);

			#GT=0/1;DP=60;AF=0.051;DLOOP:CR=1.056772;CP=10.98

			$ANNOTATION="SM=$F[-1]";
			($F[7]=~/(.*)GT=(.+?);DP=(\d+);AD=(\d+);AF=(\d+\.\d+)(.*)/ or $F[7]=~/(.*)GT=(.+?);DP=(\d+);AD=(\d+);AF=(\d+)(.*)/) or die "ERROR: $_";
			{
				$GT=$2;
				$DP=$3;
				$AD=$4;
				$AF=$5;

				$ANNOTATION.=";$1" if($1); 
				$ANNOTATION.="$6"  if($6);
				$ANNOTATION=~s/;;/;/;
			}

			$F[7]=$ANNOTATION;
			$F[8]="GT:DP:AD:AF";
			$F[9]="$GT:$DP:$AD:$AF";
		}

		print join "\t",@F;
		print "\n";
	}
	exit 0;
}

