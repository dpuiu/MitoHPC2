#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that filters unique SNV's (same ref,pos,ref,alt,sample)

        EXAMPLE:
                cat I.vcf -min min_count | $0 

~;

###############################################################################
#
# Main program
#
###############################################################################


MAIN:
{
        my %opt;
	$opt{min}=1;
        my $result = GetOptions(
                "help"          => \$opt{help},
		"min=i"		=> \$opt{min}
        );
        if(!$result)            { die "ERROR: $! "}
        if($opt{help})          { print $HELP; exit 0 }

        ##########################################################################

	my (%header,%line,%count);
	while(<>)
	{
		if(/^#/)
		{
			print unless($header{$_}); 
			$header{$_}=1;
			next;
		}

		my @F=split /\t/;

		my $SM="";
		my $AF=1;
		if(@F>8)
		{
			if($F[8] eq "SM") { $SM=$F[9]}
			elsif($F[7]=~/SM=(.+?);/ or $F[7]=~/SM=(.+)$/) { $SM=$1}
		}
	
		my $key=join "\t",(@F[0..4],$SM);
		$line{$key}=$_ unless($line{$key});		
		$count{$key}++;
	}

	foreach my $key ( keys %line )
	{
		print $line{$key} if($count{$key}>=$opt{min});
	}
}

