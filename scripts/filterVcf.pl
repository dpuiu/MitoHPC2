#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that filters a VCF file

        EXAMPLE:
                 cat I.vcf | $0 -sample s -source m -p p > O.vcf

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
	$opt{percent}="0.0";
	$opt{depth}=0;
	my %suspicious;

	my $result = GetOptions(
                "percent=s" 	=> \$opt{percent},
		"depth=i"	=> \$opt{depth},
		"sample|Run=s"	=> \$opt{sample},
		"source=s"	=> \$opt{source},
		"header=s"	=> \$opt{header},
		"suspicious=s"	=> \$opt{suspicious},
		"help"  	=>\$opt{help}
        );
        if(!$result)            { die "ERROR: $! "}
        if($opt{help})          { print $HELP; exit 0 }
	$opt{percent}=~/([-+]?([0-9]+(\.[0-9]+)?|\.[0-9]+))/ and $opt{percent}>=0 and $opt{percent}<=1  or die "ERROR:percent must be >=0 and <=1";

	############################################################################

	if($opt{header})
	{
		open(IN,$opt{header}) or die "ERROR: $!";
		print while(<IN>);
		close(IN);
			
		print "##sample=$opt{sample}\n" if(defined($opt{sample}));
		print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
	}
	
	if($opt{suspicious})
	{
                open(IN,$opt{suspicious}) or die "ERROR: $!";
                while(<IN>)
                {
                        my @F=split /\t/;
			$suspicious{$F[0]}=1 if(@F);
                }
                close(IN)
        }

	while(<>)
	{
		if(/^#/)
		{
			print unless($opt{header}) ;
			next ;
		}

		chomp;
                my @F=split /\t/;
		die "NORMALIZATION ERROR:$_\n" if($F[4]=~/,/);
		#$F[6]="." unless($F[6] eq "PASS");

		if($F[8] eq "SM")
		{
			if($F[7]=~/(.*)AF=(0\.\d+)(.*)/ or $F[7]=~/(.*)AF=(\d)(.*)/)
			{
				my $AF=$2;
				if($AF<$opt{percent})            { next; }
				elsif($AF>1-$opt{percent})	 { $AF=1; }

				$F[7]="$1AF=$AF$3";
			}
			if($F[7]=~/DP=(\d+)/ and $1<$opt{depth}) { next }

			next if($suspicious{$F[9]});
			print join "\t",@F[0..9];
			print "\n";
		}
		else
		{
	                if($F[7]!~/SM=/ and defined($opt{sample}))
        	        {
                	        $F[7]="SM=$opt{sample}";
				$F[7].=";INDEL" if($F[7]!~/INDEL/ and ($F[4] eq "*" or length($F[3]) ne length($F[4])));
        	        }

			my @F8=split /:/,$F[8];
			my @F9=split /:/,$F[9];
			my %h;
			foreach my $i (0..@F8-1)
			{
				$h{$F8[$i]}=$F9[$i];
			}
			$h{AF}=1 unless(defined($h{AF}));
			$h{AF}=$1 if($h{AF}=~/(\d\.\d+)/);

			if(!defined($h{AF}) or $h{AF}>1-$opt{percent}) 	{ $h{AF}=1; }
			elsif($h{AF}<$opt{percent})			{ next;     }

			$h{GT}="1$1" if($h{AF}==1 and $h{GT}=~/0(.+)/);

			if($h{DP}<$opt{depth}) { next }

			$F[8]="GT:DP:AF";
			$F[9]="$h{GT}:$h{DP}:$h{AF}";


			next if($F[7]=~/SM=([^;\s]+)/ and $suspicious{$1});
			print join "\t",@F[0..9];
			print "\n";
		}
	}
	exit 0;
}
