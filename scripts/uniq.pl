#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that filters up to max(1) lines with uniq values in a certain column

        EXAMPLE:
                cat I.tab  | $0 -i N -max M

~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
	my $i;
	my $max=1;
	my %opt;

        my $result = GetOptions(
                "i=i"   	=> \$i,
		"count|max=i"	=> \$max,
 		"help"  	=> \$opt{help}
        );
        if(!$result)            { die "ERROR: $! "}
        if($opt{help})          { print $HELP; exit 0 }

	##########################################################################
	my %h;
	while(<>)
	{
		chomp;
		if(defined($i))
		{
			if(/^#/ or /^$/) { print "$_\n"}
			elsif(/^@/) { print "$_\n"}
			else
			{				
				my @F=split /\t/;
				die if(!defined($F[$i]));

				$h{$F[$i]}++;
				print "$_\n" if($h{$F[$i]}<=$max);
			}
		}
		else
		{
			$h{$_}++;
                        print "$_\n" if($h{$_}<=$max);
		}
	}

	exit 0;
}

