#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that filters a SAM file ; both the reads and its mate must be aligned in the intervals given as arguments

        EXAMPLE:
                 samtools view -h I.bam  | $0 chrM:1-16569 chr1:629084-634672 chr17:22521208-22521639 | samtools view b > O.bam

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
	my $result = GetOptions(
		"help"  =>\$opt{help}
	);
	if(!$result)            { die "ERROR: $! "}
        if($opt{help})          { print $HELP; exit 0 }

	#######################################################################
	my (%h,%hb,%he);

	foreach(@ARGV)
	{
		if($_=~/(\S+):(\d+)-(\d+)/) 
		{ 
			push @{$hb{$1}},$2 ; 
			push @{$he{$1}},$3 ; 
		}
		elsif($_=/(\S+)/)
		{
			$h{$1}=1
		}
	}

	while(<>)
	{
		chomp;
		my @F=split /\t/;

		if(/^\@SQ\tSN:(\S+)/ and $hb{$1})
		{
			print "$_\n";
		}
		elsif(/^\@SQ\tSN:(\S+)\s+LN:(\d+)/ and $h{$1})
                {
			push @{$hb{$1}},1 ;
			push @{$he{$1}},$2; 
                        print "$_\n";
                }
		elsif(/^\@SQ\t/)
		{
		}
		elsif(/^\@/) 
		{
			print "$_\n";
		}
		elsif(@F>7 and $hb{$F[2]} and ($F[6] eq "=" or $hb{$F[6]}))
		{
			my $keep=0;
			foreach my $i (0..scalar(@{$hb{$F[2]}})-1)
			{
				if($hb{$F[2]}[$i]<=$F[3] and $F[3]<=$he{$F[2]}[$i]) 
				{ 
					$keep++; 
					last 
				}
			}

			$F[6]=$F[2] if($F[6] eq "=");
                        foreach	my $i (0..scalar(@{$hb{$F[6]}})-1)
                        {
                                if($hb{$F[6]}[$i]<=$F[7] and $F[7]<=$he{$F[6]}[$i]) 
				{ 
					$keep++; 
					last 
				}
                        }

			if($keep==2)
			{
				$F[9]="*"; $F[10]="*";
				print join "\t",@F; print "\n";
			}
		}
	}
	exit 0;
}
