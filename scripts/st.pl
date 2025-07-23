#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that generates summary statistics

        EXAMPLE:
                cut -f3 I.cvg  | $0 > O.cvg.stat

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
	my (@c,$s,$n);

        my $result = GetOptions(
 		"sample|Run=s"	=> \$opt{sample},
		 "help" 	=> \$opt{help}
	);
        if(!$result)            { die "ERROR: $! "}
        if($opt{help})          { print $HELP; exit 0 }

	#######################################################################

        while(<>)
        {
		chomp;
		/^\d+$/ or /^\d+\.\d+$/ or die "ERROR: $_";
		push @c,$_;
		$s+=$_;
        }

	@c=sort {$a<=>$b} @c;


	if(@c and @c>1)
	{
		$n=@c;
		print("Run\t") if(defined($opt{sample}));
		print join "\t",("count","min","q1","median","q3","max","mean"); print "\n";

		print("$opt{sample}\t") if(defined($opt{sample}));
		print join "\t",(scalar(@c),$c[0],$c[int($n/4+.5)],$c[int($n/2+.5)],$c[int($n*3/4+.5)],$c[-1],int($s*1000/$n+.5)/1000); print "\n";
	}

	exit 0;
    }

