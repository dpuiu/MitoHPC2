#!/usr/bin/env perl 
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that corrects deepvariant VCF output

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

		$F[8]=~/^GT:GQ:DP:AD:VAF:/ or die "ERROR: $_";
		$F[9]=~/^(.+?):.+?:(\d+):(\d+,\d+):(.+?):/ or die "ERROR: $_";
                my ($GT,$DP,$AD,$AF)=($1,$2,$3,$4);
		$GT="0/1" if($GT eq "./." or $GT eq "0/0");
                $F[8]="GT:DP:AD:AF";
                $F[9]="$GT:$DP:$AD:$AF";

		print join "\t",@F; print "\n";
	}
}

