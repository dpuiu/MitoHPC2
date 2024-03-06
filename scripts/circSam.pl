#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that "circularizes" SAM alignments

	EXAMPLE:
		 cat I.sam | $0 -ref_len R.fa.fai -offset N > O.sam

~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
	# define variables
	my %options;
	my %len;
	$options{offset}=0;

	# validate input parameters
	my $result = GetOptions(
		"ref_len=s"	=>	\$options{ref_len},
		"offset=i"	=>	\$options{offset},
		"help" 	        =>	\$options{help}
	);

        if(!$result)                    { die   "ERROR: $! "}
	if($options{help})              { print $HELP; exit 0 }
        if(!defined($options{ref_len})) { die   "ERROR: missing reference file" }
        
	open(IN,$options{ref_len}) or die;
        while(<IN>)
	{
		chomp;
                next if(/^#/ or /^$/);
             	next unless(/^(\S+)\s+(\d+)/);

        	$len{$1}=$2;
	}
	close(IN);
	########################################
	while(<>)
	{
		chomp;
		my @F=split /\t/;
		if(/^\@SQ\s+SN:(\S+)\s+LN:(\S+)/)
		{
			if($len{$1})
			{
				$F[2]="LN:$len{$1}";
			}
		}
		elsif(/^\@/)
                {
                }
		elsif($F[5]=~/M/)	# new
		{
			$F[3]+=$options{offset};
			$F[7]+=$options{offset};

			if($F[7])
			{
				if($F[6] eq "=" and $len{$F[2]} and $F[7]>$len{$F[2]})
				{
					$F[7]%=$len{$F[2]}
				}
				elsif($len{$F[6]} and $F[7]>$len{$F[6]})
                               	{
                                       	$F[7]%=$len{$F[6]}
                               	}
			}

			if($F[3] and $len{$F[2]} and $F[3]<=$len{$F[2]})
			{
				my ($cigar1,$cigar2)=("","");
				my ($trim1,$trim2)=(0,0);
				my $pos=$F[3]-1;

				if($F[5]=~/^(\d+)([SH])(.+)/) 
				{
					$cigar1=$1.$2;
					$trim2+=$1;
					$F[5]=$3;
				}

				while($F[5]=~/^(\d+)(\w)(.*)/)
				{
					if($pos+$1<=$len{$F[2]})
					{
						$cigar1.=$1.$2;
						$trim2+=$1 unless($2 eq "D");
					}
					elsif(!$trim1 and $2 eq "D")
					{
						$cigar1.=$1.$2;
					}
					elsif(!$trim1)
					{
						$trim1=$pos+$1-$len{$F[2]};
						my $keep1=$1-$trim1;

						if($keep1)
						{
							$trim2+=$keep1;
							$cigar1.=$keep1.$2;
						}

						$cigar2=$trim2."S".$trim1.$2;
					}
					else
					{
						$cigar2.=$1.$2;
						$trim1+=$1 unless($2 eq "D");
					}

					#$pos+=$1 unless($2 eq "D");  #June 24
					$pos+=$1 unless($2 eq "I");

					$F[5]=$3;
				}

				if($trim1)
				{
					$cigar1.=$trim1."S";
				}

				if($cigar1=~/^(.+?)(\d+)S(\d+)S$/)
				{
					$cigar1=$1.($2+$3)."S";
				}

				if($cigar2=~/M/)
				{

					print join "\t",($F[0],$F[1]+2048,$F[2],1,$F[4],$cigar2,@F[6..scalar(@F)-1]); # June 24 2020
					#print join "\t",($F[0],$F[1]+2048,$F[2],1,$F[4],$cigar2,@F[6..10],"NM:i:0",@F[12..scalar(@F)-1]);  #aug 14, changed from 0 to 1; temp
					print "\n";
				}

			        if($cigar1=~/M/)
                               	{
					$F[5]=$cigar1;
					#$F[11]="NM:i:0" ;  # temp; only on example3
                               	}
				else
				{
					@F=()
				}

			}
		}

		if(!/^@/ and $F[3]=~/^\d+$/ and $len{$F[2]} and $F[3]>$len{$F[2]})
		{
			$F[3]=$F[3]%$len{$F[2]};
		}
		
		print join "\t",@F if(@F);
		print "\n";
	}
	exit 0;
}

