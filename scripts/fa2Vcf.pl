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
	if(@ARGV<1)             { die "ERROR: insufficient number of arguments"}

	######################################################################

	print "##reference=file://$ARGV[0]\n";

	open(IN,"$ARGV[0].fai") or die "ERROR:$!";
	while(<IN>)
	{
		chomp;
                my @F=split /\t/;
		print "##contig=<ID=$F[0],length=$F[1]>\n";
	}
	print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";

	exit 0;
}

