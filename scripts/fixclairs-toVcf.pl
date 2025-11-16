#!/usr/bin/env perl 
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that corrects mutect2 VCF output

        EXAMPLE:
                 cat I.vcf | $0 > O.vcf

~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
        my %opt;
	my (%h,%AF);
        my $result = GetOptions(
		"file=s"         => \$opt{file},
                "help"           => \$opt{help}
        );
	if(!$result)            { die "ERROR: $! "}
        if($opt{help})          { print $HELP; exit 0 }

	############################################################

	while(<>)
	{

		if(/^#/)
		{
			print;
			next
		}
		chomp;
		my @F=split /\t/;
		next if($F[3]=~/N/);

		$F[8]=~/^GT:GQ:DP:AF:AD/ or die "ERROR: $_";
		$F[9]=~/^(.+?):.+?:(.+?):(.+?):(.+?):/ or die "ERROR: $_";
                my ($GT,$DP,$AF,$AD)=($1,$2,$3,$4);
                $F[8]="GT:DP:AD:AF";
                $F[9]="$GT:$DP:$AD:$AF";

		print join "\t",@F; print "\n";
	}
}

