#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that merges 2 POS files

        EXAMPLE:
                cat I.pos | $0  > O.pos

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

        #########################################

	#0     1    2   3    4    5	  						   6       7	   8	   9	   10	   11	  
	#chr   pos  ID  ref  alt  filters                                                  AC_hom  AC_het  AF_hom  AF_het  AN      max_ARF
	#chrM  12   .   T    C    DLOOP;MLC_score=0.56635886;MLC_consq=intergenic_variant  6       0       3e-05   0       199955  0

	my $AN=0;
	while(<>)
	{
		if(/^chr\t/) { print if($.==1); }
		else
		{
			chomp;
			my @F=split /\t/;
			my $key= join "\t",@F[0,1,3,4];
			$h{$key}{info}=join "\t",@F[0..5];
			$h{$key}{AC_hom}+=$F[6];
			$h{$key}{AC_het}+=$F[7];
			$h{$key}{AN}+=$F[10];
			$AN=$h{$key}{AN} if($AN<$h{$key}{AN});
			$h{$key}{max_ARF}=$F[11] if(!$h{$key}{max_ARF} or $h{$key}{max_ARF}<$F[11]);
		}
	}
	
	foreach my $key (keys %h)
	{
		print join "\t",($h{$key}{info},$h{$key}{AC_hom},$h{$key}{AC_het},$h{$key}{AC_hom}/$AN,$h{$key}{AC_het}/$AN,$AN,$h{$key}{max_ARF});
		print "\n";
	}
		
	exit 0;
}

