#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that filters unique SNV's (same ref,pos,ref,alt,sample)

        EXAMPLE:
                cat I.vcf  | $0 

~;

###############################################################################
#
# Main program
#
###############################################################################


MAIN:
{
        my %opt;
        my $result = GetOptions(
                "help"          => \$opt{help}
        );
        if(!$result)            { die "ERROR: $! "}
        if($opt{help})          { print $HELP; exit 0 }

        ##########################################################################

	my (%max,%line);
	while(<>)
	{
		if(/^#/)
		{
			print; next;
		}

		my @F=split /\t/;

		my $SM="";
		my $AF=1;
		if(@F>8)
		{
			if($F[8] eq "SM") { $SM=$F[9]}
			elsif($F[7]=~/SM=(.+?);/ or $F[7]=~/SM=(.+)$/) { $SM=$1}

			$AF=$1 if($F[7]=~/AF=(\S+?);/ or $F[7]=~/AF=(\S+)$/);
		}
	
		my $key=join "\t",(@F[0..4],$SM);

		if(!$max{$key} or $AF>$max{$key})
		{
			$max{$key}=$AF;
			$line{$key}=$_;
		}
	}

	foreach my $key ( keys %max )
	{
		print $line{$key};
	}
}

