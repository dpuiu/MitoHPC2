#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that converts a concatenated VCF file (one SNV/line) to a merged VCF file (one column/sample); 

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

	#######################################################

	my ($pkey,$sample,@keys,%keys,%GT_DP_AF,$header);
	my ($AC,$AN)=(0,0);

	while(<>)
	{
		if(/^##/)
		{
			next if(/^##FILTER/);
			print;
		}
		elsif(/^#/)
		{

			print "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n";
			print "##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Allele count in genotypes\">\n";
		}
		else
		{
			#0	1	2	3	4	
			#chrM	42	.	T	TC		

			if(!$header) { print join "\t",("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",@samples); print "\n"; $header=1}

			chomp;
			my @F=split /\t/;
			@F[2,5,6]=(".",".",".");
	
			if($F[-1]=~/^(.+):(.+):1$/)        { $F[-1]="1:$2:1"    }
			elsif($F[-1]=~/^(.+):(.+):(.+)$/)  { $F[-1]="0/1:$2:$3" }
	

			my $sample;
			if($F[7]=~/SM=(.+?);(.+)/) { ($sample,$F[7])=($1,"$2;");}
			elsif($F[7]=~/SM=(.+)/)    { ($sample,$F[7])=($1,"");}

			defined($samples{$sample}) or die "ERROR: sample not defined in $_";
			my $key=join "\t",@F[0..7];
			
			if($pkey and $key ne $pkey)
			{
 				my @P=($pkey."AC=$AC;AN=$AN","GT:DP:AF");

                		foreach my  $sample (@samples)
		                {
		                        if($GT_DP_AF{$sample}) { push @P,$GT_DP_AF{$sample} ; }
                	       	 	else                   { push @P,".:.:." }
	                	}

	        	        print join "\t",@P;  print "\n";

				%GT_DP_AF=();
				($AC,$AN)=(0,0);
				
			}
			
			$pkey=$key;

			$GT_DP_AF{$sample}=$F[-1];			
			$AC++;
                        $AN++;
                        $AN++ if($GT_DP_AF{$sample}=~/\|/ or $GT_DP_AF{$sample}=~/\//)

		}
	}

	#end
	if($pkey)
        {
		my @P=($pkey."AC=$AC;AN=$AN","GT:DP:AF");

		foreach my  $sample (@samples)
                {
			if($GT_DP_AF{$sample}) { push @P,$GT_DP_AF{$sample} ; }
                        else                   { push @P,".:.:." }
                }

               print join "\t",@P;  print "\n";
	}

	exit 0;
}
