#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that converts a merged VCF file(one column/sample) into a concatenated VCF file (one SNV/line)

        EXAMPLE:
                cat I.vcf | $0 -in I.txt [ -suspicious S.txt ]  > O.merge.vcf

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
	my %suspicious;

	my $result = GetOptions(
		"in=s"		=> \$opt{in},
                "suspicious=s"  => \$opt{suspicious},
		"help" 		=> \$opt{help}
	);
        if(!$result)    { die "ERROR: $! "}
        if($opt{help})  { print $HELP; exit 0 }

	#######################################################

        if($opt{suspicious})
        {
                open(IN,$opt{suspicious}) or die "ERROR: $!";
                while(<IN>)
                {
			chomp;
                        my @F=split /\t/;
                        $suspicious{$F[0]}=1 if(@F);
                }
                close(IN)
        }

        my (@samples,%samples);
	if($opt{in})
	{
		open(IN,$opt{in}) or die "ERROR: $!";
        	while(<IN>)
	        {
			chomp;
			next if(/^#/ or /^$/);
			my @F=split /\t/;
        	        next if($suspicious{$F[0]});
                	push @samples,$F[0];
			$samples{$F[0]}=1;
        	}
		close(IN);
	}

	#######################################################

	while(<>)
	{
		if(/^##/)
		{
			print;
		}
		elsif(/^#/)
		{
			print join "\t",("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE"); 
			print "\n"; 

			if(!@samples)
			{
				my @F=split;
				@samples=@F[9..@F-1];
			}
		}
		else
		{		
			chomp;
			my @F=split /\t/;
	
			foreach my $i (9..@F-1)
			{
				if($F[$i]=~/\d/)
				{
					print join "\t",(@F[0..6],"$F[7];SM=".$samples[$i-9],$F[8],$F[$i]);
					print "\n";
				}			
			}
		}
	}

	exit 0;
}
