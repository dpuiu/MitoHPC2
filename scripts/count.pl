#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that counts the values in a specific column (default 0)

        EXAMPLE:
                 cat I.tab | $0 -i i -j j  > O.count

~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
	my $i=0;
	my $j;
	my $min;
	my $max;
	my %count;
	my $round;
	my $help;

        my $result = GetOptions(
		"i=i"		=> \$i,
		"j=i"		=> \$j,
		"min=i"		=> \$min,
		"max=i"		=> \$max,
		"round=i"	=> \$round,
		"help"		=> \$help
	);
	if(!$result)    { die "ERROR: $! "}
        if($help)  	{ print $HELP; exit 0 }

	######################################################
	while(<>)
	{
		chomp;
                next if(/^$/ or /^#/) ;

		my @F=split /\t/;
		#my @F=split;

		next if(@F<=$i);
		next if(defined($j) and @F<=$j);

		$F[$i]=int($F[$i]*$round+.5)/$round if(defined($round));

		if(defined($j)) { $count{$F[$i]}+=$F[$j] }
		else		{ $count{$F[$i]}++;      }
	}

	##########################################################

	foreach my $key (sort {$count{$b}<=>$count{$a}} keys %count)
	{
		next if(defined($min) and $count{$key}<$min);
		next if(defined($max) and $count{$key}>$max);

                print $key,"\t", $count{$key};
                print "\n";
	}

	exit 0;
}

