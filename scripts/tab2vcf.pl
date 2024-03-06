#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that converts a TAB to  a VCF file

	EXAMPLE:
		$0 I.tab > O.vcf
		cat I.tab | $0 > O.vcf  

	IN:
		#0    1	     2             3    4    5                 6                 7           
		Chr   Start  dbSNP_155_id  Ref  Alt  HelixMTdb_AC_hom  HelixMTdb_AF_hom  HelixMTdb_AC_het  HelixMTdb_AF_het  HelixMTdb_mean_ARF  HelixMTdb_max_ARF  Gnomad31_filter  Gnomad31_AC_hom  Gnomad31_AC_het  Gnomad31_AF_hom  Gnomad31_AF_het  Gnomad31_AN  MITOMAP_Disease_Clinical_info  MITOMAP_Disease_Status
		chrM  3307   .             A    G    1.0               5.1024836e-06     2.0               1.0204967e-05     0.27253             0.31429            .                .                .                .                .                .            .                              .
		chrM  3307   rs1603218882  A    T    .                 .                 .                 .                 .                   .                  .                .                .                .                .                .            .                              .

	OUT
		#0      1       2       3       4       5       6       7 
		#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
		chrM    3307	.	A	G	.	.	HelixMTdb_AC_hom=1.0;HelixMTdb_AC_het=2.0;



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
	my $result = GetOptions();
        if(!$result)            	{ die "ERROR: $! "}
        if($opt{help})		{ print $HELP; exit 0 }

	###############################################################################

        print "##fileformat=VCFv4.2\n";
	print "##reference=chrM.fa\n";
	print "##contig=<ID=chrM,length=16569>\n";

	my (@H,@T,$O);
	while(<>)
	{
		chomp;
		s/;/,/g;
		my @F=split /\t/;

		if($.==1)
		{			
			@H=@F;

			foreach  (6..@H-1)
			{
				$T[$_]="String";
			}
		}
		else
		{
			my @O=(@F[0..4],".",".");
			$O[7]="";
			foreach  (6..@F-1)
			{
				next if($F[$_] eq ".");

				if($F[$_]=~/^[-+]?[0-9]+$/)                     { $T[$_]="Integer" }
				elsif($F[$_]=~/^[-+]?\d*[.]?\d*[eE]?[-+]?\d*$/) { $T[$_]="Float"   }

				if($T[$_] eq "String") { $O[7].="$H[$_]=\"$F[$_]\";"; }
				else  		       { $O[7].="$H[$_]=$F[$_];";     }
			}
			if($O[7]) 
			{ 
				chop $O[7];
				my $OO=join "\t",@O;
				$OO.="\n";
				$O.=$OO;
			}
		}
	}
                       
	foreach  (6..@H-1)
	{
		print "##INFO=<ID=$H[$_],Number=1,Type=$T[$_],Description=\"$H[$_]\">\n";
        }

        print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
	print $O;

	exit 0;
}
