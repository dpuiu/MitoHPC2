#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that counts SNVs
        H:Homoplasmies
        h:heteroplasmies
        HS:Homoplasmies(no INDELs)
        hS:heteroplasmies(no INDELs)
        HI:Homoplasmies(INDELs)
        hI:heteroplasmies(INDELs)

        EXAMPLE:
                cat I.vcf | $0 > O.count

~;

###############################################################################
#
# Main program
#
################################################################################

MAIN:
{
	# define variables
	my %opt;

	my $result = GetOptions(
		"total=i"	=> \$opt{total},
		"help"  	=> \$opt{help}
	);
        if(!$result)             { die "ERROR: $! "}
        if($opt{help})           { print $HELP; exit 0 }

	#########################################################################

	my %pos;
	while(<>)
	{
		next if(/^#/);
		chomp;
		my @F=split /\t/;
		my $t=(/AF=0/ or /:0\.\d+$/)?"h":"H";

		$pos{"$F[0]\t$F[1]"}{$t}++;
		$t.=(/INDEL/)?"I":"S";

		$pos{"$F[0]\t$F[1]"}{$t}++;
	}

	print join "\t",("#chr","pos","H","h","HS","hS","HI","hI"); print "\n";
	foreach my $k (keys %pos)
	{
		foreach my $t ("H","h","HS","hS","HI","hI")
		{
			$pos{$k}{$t}=0 unless($pos{$k}{$t});
		}

                foreach my $t ("H","h","HS","hS","HI","hI")
                {
			
                        $pos{$k}{$t}=int(10000*$pos{$k}{$t}/$opt{total}+.5)/100 if($opt{total});
                }

		print join "\t",($k,$pos{$k}{"H"},$pos{$k}{"h"}, $pos{$k}{"HS"},$pos{$k}{"hS"},$pos{$k}{"HI"},$pos{$k}{"hI"});
		print "\n";
	}
	exit 0;
}

