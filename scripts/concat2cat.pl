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
        if($opt{help})	{ print $HELP; exit 0 }

	my %GT;
	while(<>)
	{
		if(/^#/)
		{
                        s/FORMAT/INFO/ if(/^##FORMAT=<ID=DP,/ or /##FORMAT=<ID=AF,/);
                        s/INFO/FORMAT/ if(/^##INFO=<ID=SM,/); 

			print "##INFO=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" if(/^##INFO=<ID=DP,/);
			print;
			next;
		}

		chomp;
		my @F=split /\t/;
		if($F[8] ne "SM")
		{
		
			my ($SM,$ANNOTATION);
			if($F[7]=~/(\S*)SM=(\S+?);(\S+)/) { ($SM,$ANNOTATION)=($2,"$1$3;") }
			elsif($F[7]=~/(\S*)SM=(\S+)/)     { ($SM,$ANNOTATION)=($2,"$1") }
			else                              { die "ERROR: $_";           }

			my ($GT,$DP,$AF)=split /:/,$F[-1];

			$F[7]=$ANNOTATION."GT=$GT;DP=$DP;AF=$AF";
			$F[8]="SM";
			$F[9]="$SM";
		}

		print join "\t",@F;
		print "\n";
	}
	exit 0;
}

