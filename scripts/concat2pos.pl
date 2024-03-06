#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that converts a concatenated VCF file (one SNV/line) to a POS file ();

        EXAMPLE:
                cat I.vcf | $0 -in I.txt [ -suspicious S.txt]  > O.pos

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
	my %h;
	my %suspicious;

	my $result = GetOptions(
		"in=s"		=> \$opt{in},
		"suspicious=s"  => \$opt{suspicious},
		"help"  	=> \$opt{help}
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


        my $AN=0;
	open(IN,$opt{in}) or die "ERROR: $!";
        while(<IN>)
        {
		chomp;
		next if(/^#/ or /^$/);
		my @F=split /\t/;
		next if($suspicious{$F[0]});
		$AN++;
        }
	close(IN);

	#######################################################
	my ($key,$sample,@keys,%keys,%GT_DP_AF);
	while(<>)
	{
		if(/^#/)
		{
		}
		else
		{
			chomp;
			my @F=split /\t/;
			my $key="$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]";

			if($F[7]=~/^SM=.+?;(.+)/)
			{
				$h{$key}{filter}=$1;
			}
			else
			{
				$h{$key}{filter}="."
			}

			if(/.+:1$/)
			{
				$h{$key}{AC_hom}++;
			}
			elsif(/.+:(0\.\d+)$/)
                        {
                                $h{$key}{AC_het}++;
				my $AF=$1;
				$h{$key}{max_observed_heteroplasmy}=$AF if(!$h{$key}{max_observed_heteroplasmy} or $AF>$h{$key}{max_observed_heteroplasmy});
			}
			elsif(/.+:0$/)
			{
				next
			}
			else
			{
				die "ERROR: $_\n"
			}

  			if(!$keys{$key})
                       	{
                               	$keys{$key}=1;
                               	push @keys,$key;
                        }

		}
	}

	print "chr\tpos\tID\tref\talt\tfilters\tAC_hom\tAC_het\tAF_hom\tAF_het\tAN\tmax_observed_heteroplasmy\n";
	foreach my $key (@keys)
	{
		my @F;

		foreach ("AC_hom","AC_het","max_observed_heteroplasmy")
		{
			$h{$key}{$_}=0 unless(defined($h{$key}{$_}));
		}

                $h{$key}{AF_hom}=int(200000*$h{$key}{AC_hom}/$AN+.5)/200000;
		$h{$key}{AF_het}=int(200000*$h{$key}{AC_het}/$AN+.5)/200000;


		push @F,($key,$h{$key}{filter},$h{$key}{AC_hom},$h{$key}{AC_het},$h{$key}{AF_hom},$h{$key}{AF_het},$AN,$h{$key}{max_observed_heteroplasmy});
		print join "\t",@F; 
		print "\n";
	}
	exit 0;
}

