#!/usr/bin/perl -w
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that identifies converts SAM to BED

	EXAMPLE:
		cat I.sam | $0 > O.bed

~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
	# define variables
	my %options;
	my %len;

	# validate input parameters
	my $result = GetOptions(
                "help"  => \$options{help}
		);
        if(!$result)             { die "ERROR: $! "}
        if($options{help})       { print $HELP; exit 0 }

	####################################

	while(<>)
	{
		next if(/^\@/);
		chomp;
		my @F=split /\t/;
                next unless(@F);		
		$F[3]--;

		while($F[5] and $F[5]=~/^(\d+)(.)(.*)/)
		{				
			if($2 eq "M" or $2 eq "=" or $2 eq "X" or $2 eq "D")
			{	                       		
				print join "\t",($F[2],$F[3],$F[3]+$1,"$1$2",$1,"+"); print "\n";
				$F[3]+=$1;
			}
                        elsif($2 eq "I")
                        {
                                        
				print join "\t",($F[2],$F[3],$F[3],"$1$2",0,"+"); print "\n";
                        }
                        $F[5]=$3;
		}
	}
	exit 0;
}

