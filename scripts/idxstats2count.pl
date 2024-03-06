#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;


my $HELP = qq~
Program that counts the total and the chrM alignments from a samtools idxstats file

        EXAMPLE:
                cat I.idxstats | $0 -sample s -chrM chrM > O.count

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
	my %count;
	$opt{sample}=".";
	$opt{chrM}="chrM";

	my $result = GetOptions(
		"sample|Run=s"	=> \$opt{sample},
		"chrM=s"	=> \$opt{chrM},
		"help"   	=> \$opt{help}
	);

	if(!$result)             { die "ERROR: $! "}
        if($opt{help})           { print $HELP; exit 0 }
	
	#######################################################################

	while(<>)
	{
		chomp;
		my @F=split /\t/;
		$count{all}+=$F[2]+$F[3];
		$count{mapped}+=$F[2];
		$count{chrM}=$F[2] if($F[0] eq $opt{chrM}); 
	}

	print join "\t",("Run","all_reads","mapped_reads","MT_reads"); print "\n";
	print join "\t",($opt{sample},$count{all},$count{mapped},$count{chrM}); print "\n";

	exit 0;
}

