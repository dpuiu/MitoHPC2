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

        # validate input parameters
        my $result = GetOptions(
		"sm"      =>      \$opt{sm},
		"help"	  =>	  \$opt{help},
		"pos"   =>	\$opt{pos},
                "snp"    =>     \$opt{snp},
                "ins"    =>     \$opt{ins},
                "del"    =>     \$opt{del},
                "het"    =>     \$opt{het},
                "hom"    =>     \$opt{hom},

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

                if($opt{snp})   { next if(length($F[3]) ne length($F[4])) }
                if($opt{ins} )  { next if(length($F[3]) ge length($F[4])) }
                if($opt{del})   { next if(length($F[3]) le length($F[4])) }
                next if($opt{het} and !(/\tAF=0/ or /;AF=0/ or /:0\.\d+$/));
                next if($opt{hom} and (/\tAF=0/ or /;AF=0/ or /:0\.\d+$/));


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

		my $key;
                if($opt{pos}) { $key="$F[0] $F[1] $SM"}
                else          { $key="$F[0] $F[1] $F[3] $F[4] $SM" }
                $h{$key}=$AF;

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

 	        next if($F[1]==3105);	# new 
                next if($F[1]==3106);	# new
                next if($F[3]=~/N/);	# new

                if($opt{snp})   { next if(length($F[3]) ne length($F[4])) }
                if($opt{ins} )  { next if(length($F[3]) ge length($F[4])) }
                if($opt{del})   { next if(length($F[3]) le length($F[4])) }
                next if($opt{het} and !(/\tAF=0/ or /;AF=0/ or /:0\.\d+$/));
                next if($opt{hom} and (/\tAF=0/ or /;AF=0/ or /:0\.\d+$/));

               	my $SM="";

		if($opt{sm})
		{
	                if(@F>8 and $F[8] eq "SM")                                   { $SM=$F[9] }
        	        elsif(@F>7 and ($F[7]=~/SM=(\S+?);/ or $F[7]=~/SM=(\S+?)$/)) { $SM=$1 }
			else                                                         { die "ERROR: $_" }
		}

		($F[3],$F[4])=(uc($F[3]),uc($F[4]));

		my $key;	# 2025/11/14

                if($opt{pos}) { $key="$F[0] $F[1] $SM"}
                else          { $key="$F[0] $F[1] $F[3] $F[4] $SM" }

	        print "$_\n" if $h{$key};
        }

	exit 0;
}
