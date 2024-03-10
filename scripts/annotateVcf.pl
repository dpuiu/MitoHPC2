#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that annotates a VCF file using the information (INFO field) from another VCF file

        EXAMPLE:
                $0 I.vcf A.vcf > O.vcf

~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
	my %options;
        my (%info,%h);

        # validate input parameters
        my $result = GetOptions(
		"help"	=>	\$options{help}
		);

        if(!$result)            { die "ERROR: $! "}
	if($options{help})      { print $HELP; exit 0 }
        if(@ARGV<2)             { die "ERROR: insufficient number of arguments"}

        #########################################

        open(IN,"zcat -f $ARGV[1] |") or die("ERROR: Cannot open input file".$!) ;
        while(<IN>)
        {
		##INFO=<ID=NUMT,Number=0,Type=Flag,Description="NUMT">
		#0	1	2	3	4	5	6	7
		#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
		#RSRS	15	.	C	T	.	.	NUMT

                if(/^##INFO/) { $info{$_}=""}
		elsif(/^#/)   {}
		else
		{
			chomp;
	                my @F=split /\t/;
			my $key="$F[0] $F[1] $F[3] $F[4]";
                	$h{$key}=$F[7];
		}
        }
	close(IN);
        last unless(%h);

        #########################################

        open(IN,$ARGV[0]) or die("ERROR: Cannot open input file".$!) ;

        while(<IN>)
        {
                if(/^#INFO/)
		{
			print $_;
			delete($info{$_});
		}
		elsif(/^#CHROM/)
		{
			print %info;
			print;
		}
		elsif(/^#/)
		{
			print;
		}
		else
		{	
			chomp;
	                my @F=split /\t/;

			my $key="$F[0] $F[1] $F[3] $F[4]";
			if($h{$key})
			{
				if(!$F[7]) 		 { $F[7]="$h{$key}" }
				elsif($F[7]!~/$h{$key}/) { $F[7].=";$h{$key}"; }
			}
			print join "\t",@F;
			print "\n";
		}
        }

	exit 0;
}


