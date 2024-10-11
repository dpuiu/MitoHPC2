#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that generates a VCF header corresponding to a FASTA file

        EXAMPLE:
		samtools faidx I.fa
                test -s I.fa.fai
                $0 I.fa > O.vcf

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
		"help"  =>\$opt{help}
	);
	if(!$result)            { die "ERROR: $! "}
        if($opt{help})          { print $HELP; exit 0 }

	######################################################################

	while(<>)
	{
		chomp;
                my @F=split /\t/;

		next if($.==1 or $F[2]=~/N/i or $F[3]=~/N/i); 
		my $DP=$F[4]+$F[5];
		$F[6]=~/(.+)%/;
		my $AF=$1/100; 

		@F=($F[0],$F[1],".",$F[2],$F[-1],".",".",".","GT:DP:AD:AF","0/1:$DP:$F[4],$F[5]:$AF");

		if($F[4]=~/\+(.+)/) 
		{
			($F[4],$F[7])=("$F[3]$1","INDEL")
		} 
		elsif($F[4]=~/\-(.+)/) 
		{
			($F[3],$F[4],$F[7])=("$F[3]$1",$F[3],"INDEL")
		} 
		print join "\t",@F; print "\n";
	}

	exit 0;
}

