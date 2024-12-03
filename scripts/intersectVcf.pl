#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that prints all the SNV's present in both VCF files given as arguments

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
		"sm"      =>      \$opt{sm},
		"min1=s"  =>	  \$opt{m1},
		"Max1=s"  =>      \$opt{M1},
               	"min2=s"  =>      \$opt{m2},
               	"Max2=s"  =>      \$opt{M2},
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
		($F[3],$F[4])=(uc($F[3]),uc($F[4]));

		my $SM="";

                if($opt{sm})
		{
			die unless(@F==10);
			if($F[8] eq "SM")                                 { $SM=$F[9] }
			elsif($F[7]=~/SM=(\S+?);/ or $F[7]=~/SM=(\S+?)$/) { $SM=$1 }
			else						  { die "ERROR: $_" }
		}

		if($opt{m2} or $opt{M2})
		{
			my $AF=1;
			$AF=$1 if(/AF=(0\.\d+)/ or /.+:(1)$/ or /.+:(0\.\d+)$/);

			next if(defined($opt{m2}) and $AF<$opt{m2});
			next if(defined($opt{M2}) and $AF>$opt{M2});
		}

                $h{"$F[0] $F[1] $F[3] $F[4] $SM"}=1;
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
		($F[3],$F[4])=(uc($F[3]),uc($F[4]));
               	my $SM="";

		if($opt{sm})
		{
			die "ERROR:$_\n" unless(@F==10);
			if($F[8] eq "SM")                                 { $SM=$F[9] }
                	elsif($F[7]=~/SM=(\S+?);/ or $F[7]=~/SM=(\S+?)$/) { $SM=$1 }
			else                                              { die "ERROR: $_" }
		}

		if($opt{m1} or $opt{M1})
                {
			my $AF=1;
        	        $AF=$1 if(/AF=(0\.\d+)/ or /.+:(1)$/ or /.+:(0\.\d+)$/);

               		next if(defined($opt{m1}) and $AF<$opt{m1});
	                next if(defined($opt{M1}) and $AF>$opt{M1});
		}

                print  "$_\n" if($h{"$F[0] $F[1] $F[3] $F[4] $SM"});
        }

	exit 0;
}

