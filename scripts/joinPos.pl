#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that joins 2 VCF files

        EXAMPLE:
                $0 I.vcf J.vcf > O.vcf

~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
        my (%opt,%h);

        # validate input parameters
        my $result = GetOptions(
		"help"    =>	  \$opt{help}
	);

	if(!$result)            { die "ERROR: $! "}
        if($opt{help})          { print $HELP; exit 0 }
        if(@ARGV<2)             { die "ERROR: insufficient number of arguments"}

        #########################################

        open(IN,$ARGV[1]) or die("ERROR: Cannot open input file".$!) ;
        while(<IN>)
        {
                next if(/^$/ or /^#/);
		chomp;
                my @F=split /\t/;
		
                $h{"$F[0] $F[1] $F[3] $F[4]"}=$_;
        }
	close(IN);
        last unless(%h);

        #########################################

        open(IN,$ARGV[0]) or die("ERROR: Cannot open input file".$!) ;

        while(<IN>)
        {
                if(/^$/ or /^#/) { print; next}

		chomp; 
                my @F=split /\t/;

		print $_,"\t";
		if($h{"$F[0] $F[1] $F[3] $F[4]"})
		{
			print $h{"$F[0] $F[1] $F[3] $F[4]"};
			delete $h{"$F[0] $F[1] $F[3] $F[4]"}
		}
		else
		{
			print "\t."x12;
		}
		print "\n";
        }

	foreach(keys %h)
	{
		print ".\t"x12;
		print $h{$_};
		print "\n";
	}

	exit 0;
}

