#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that computes the mtDNA-CN in each sample using aligned read counts genarated 
  by processing the idxstats file

        EXAMPLE:
                cat I.count | $0 > O.count

                #0              1          	2          	3	
		#Run		all_reads	mapped_reads	MT_reads
		#HG00438	760401407	759218290	1915191	
		#HG00621	737055935	735795660	1506529		
		
		#=>

		#0              1               2               3		4
		#Run            all_reads       mapped_reads    MT_reads        mtDNA-CN
                #HG00438        760401407       759218290       1915191         923
                #HG00621        737055935       735795660       1506529         749

~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
	# define variables
	my (%opt,%len);
	#$len{ref}=   3217346917;	#hs38DH size (3355 sequences) 
	#$len{ref}=   3099922541;	#hs38DH size (195 sequences)
	#$len{female}=3160119502;	#hs38DH  based estimates
	#$len{male}=  3110712762;	#hs38DH  based estimates

	$len{ref}=3031865587;		#CHM13 v1.1.	94.23% of prev estimate 
	$len{female}=3054815472;	#CHM13 v1.1     96.66% of prev estimate 
	$len{male}=3008915703;		#CHM13 v1.1;  	96.72% of prev estimate 

	$opt{chrM}=16569;

	my $result = GetOptions(
		"ref=i"	 =>	\$opt{ref},
		"female" =>   	\$opt{female},
		"male"	 =>     \$opt{male},
		"help"   => 	\$opt{help}
        );

	if(!$result)             { die "ERROR: $! "}
        if($opt{help})           { print $HELP; exit 0 }

	if($opt{female})   { $opt{ref}=$len{female} }
	elsif($opt{male})  { $opt{ref}=$len{male}   }
	elsif(!$opt{ref})  { $opt{ref}=$len{ref}    }

	#########################################################################
	while(<>)
	{
		chomp;
		my @F=split /\t/;
		if(@F>=4)
		{
			if($.==1 and (/^Run/ or /^sample/))
			{
				push @F,"mtDNA-CN";
			}
			else
			{
				my $M=($F[3]*2*$opt{ref})/($F[2]*$opt{chrM});
				#$M=int($M*100+.5)/100;
				$M=int($M+.5);
				$F[4]=$M ;
			}
		}
		print join "\t",@F;  
		print "\n";
	}
	exit 0;
}


