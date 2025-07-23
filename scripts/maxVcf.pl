#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that computes the alleles with AF > referece allele

        EXAMPLE:
                 cat I.vcf | $0 > M.vcf
~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{

        my %options;
        my $result = GetOptions(
                "help"  =>      \$options{help}
                );

        if(!$result)            { die "ERROR: $! "}
        if($options{help})      { print $HELP; exit 0 }

        #########################################

	my (%max,%line,%sum,%af,%line2,%af2,%line3);
	while(<>)
	{
		if(/^#/) { print; next; }

		chomp;
		my @F=split /\t/;
		my $POS=$F[1];
		($F[3],$F[4])=(uc($F[3]),uc($F[4]));

		my $AF=1;
		$AF=$1 if($F[9]=~/.+:(\S+)/);
                if(!$max{$POS} or $AF>$max{$POS})
                {
                        $max{$POS}=$AF;
                        $line{$POS}=$_;
			$af{$POS}=$AF;
                }

		$sum{$POS}+=$AF;
	}


	foreach my $POS ( keys %line )
	{
		if($max{$POS}>1-$sum{$POS})
		{
			$line2{$POS}=$line{$POS};
			$af2{$POS}=$af{$POS};
		}
	}

	my $POSMAX=0;
	foreach my $POS ( sort {$a <=> $b} keys %line2 )
        {
		 my @P=split /\t/,$line2{$POS};
		 my $DIFF=abs(length($P[3])-length($P[4]));

		 #if($POS>=$POSMAX and !($P[6]=~/multiallelic/ and $DIFF))
                 if($POS>=$POSMAX)
		 {
			print $line2{$POS},"\n";
		 }			

		 my $POSPLUS=$POS+$DIFF;
		 $POSMAX=$POSPLUS if($POSMAX<$POSPLUS);
	}
}
