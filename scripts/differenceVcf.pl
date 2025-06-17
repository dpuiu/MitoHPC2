#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that prints all the SNV's present in one VCF file but missing from another

        EXAMPLE:
                $0 I.vcf J.vcf [ -sm -af ] > O.vcf

~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
        my (%opt,%h);
	$opt{af};	

        # validate input parameters
        my $result = GetOptions(
		"sm"      =>      \$opt{sm},
               	"af=s" 	  =>      \$opt{af},
		"help"	  =>	  \$opt{help}
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
                die "ERROR $_" if(@F<5);

		my $SM="";

		if($opt{sm})
		{
			if(@F>8 and $F[8] eq "SM")                                   { $SM=$F[9] }
			elsif(@F>7 and ($F[7]=~/SM=(\S+?);/ or $F[7]=~/SM=(\S+?)$/)) { $SM=$1 }
			else						             { die "ERROR: $_" }
		}

		my $AF=1;
		$AF=$1 if(/AF=(0\.\d+)/ or /.+:(1)$/ or /.+:(0\.\d+)/);
		($F[3],$F[4])=(uc($F[3]),uc($F[4]));
                $h{"$F[0] $F[1] $F[3] $F[4] $SM"}=$AF;
        }
	close(IN);
        #last unless(%h);

        #########################################

        open(IN,$ARGV[0]) or die("ERROR: Cannot open input file".$!) ;

        while(<IN>)
        {
                if(/^$/ or /^#/) { print; next}

		chomp;
                my @F=split /\t/;
                die "ERROR $_" if(@F<5);
               	my $SM="";

		if($opt{sm})
		{
	                if(@F>8 and $F[8] eq "SM")                                   { $SM=$F[9] }
        	        elsif(@F>7 and ($F[7]=~/SM=(\S+?);/ or $F[7]=~/SM=(\S+?)$/)) { $SM=$1 }
			else                                                         { die "ERROR: $_" }
		}

		($F[3],$F[4])=(uc($F[3]),uc($F[4]));
		if($h{"$F[0] $F[1] $F[3] $F[4] $SM"})
		{
			if($opt{sm} and $opt{af})
			{
				if(/(.+;AF=)(0\.\d+)(.+)/ or /(.+;AF=)(1)(.+)/)
				{
					my $AF=$2-$h{"$F[0] $F[1] $F[3] $F[4] $SM"};
					print "$1$AF$3\n" if($AF>$opt{af});
				}	
				elsif(/(.+:)(\d.*)$/)
				{
 				 	 my $AF=$2-$h{"$F[0] $F[1] $F[3] $F[4] $SM"};
					 print "$1$AF\n" if($AF>$opt{af});
				}
			}
		}
		else
		{
	                print "$_\n" unless $h{"$F[0] $F[1] $F[3] $F[4] $SM"};
		}
        }

	exit 0;
}
