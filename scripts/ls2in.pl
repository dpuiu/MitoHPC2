#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that parses the find/ls comamnd output, identifies all the BAM/CRAM files
  and generetes the input config files required by MitoHPC

        EXAMPLE:
                find ADIR -type f -name "*.bam" | $0 -out ODIR > in.txt

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
	$opt{out}="out";

	my $result = GetOptions(
		"out=s"	=>	\$opt{out},
		 "help"  =>     \$opt{help}
	);

        if(!$result)            { die "ERROR: $! "}
        if($opt{help})          { print $HELP; exit 0 }

	###################################

	while(<>)
	{
		chomp;
		next unless(/\.bam$/ or /\.cram$/);
		my @F=split /\t/;

		print "$1\t$F[-1]\t$opt{out}/$1/$1\n" if($F[-1]=~/.+\/(\S+)\./ or $F[-1]=~/(\S+)\./);
		#print "$1\t$F[-1]\t$opt{out}/$1/$1\n" if($F[-1]=~/.+\/(\S+?)\./ or $F[-1]=~/(\S+?)\./);
	}
	exit 0;
}

