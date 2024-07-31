#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that recalculates SNV positions using the max VCF and reference FASTA

        EXAMPLE:

		cat I.vcf | $0 -mfile M.vcf -rfile R.fa  > O.vcf

~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
	my %opt;
	$opt{ref}="chrM";
	$opt{rlen}=16569;

	my $result = GetOptions(
		"ref=s"		=> \$opt{ref},		#chrM
		"rlen=i"        => \$opt{rlen},		#16569
 	        "mfile=s" 	=> \$opt{mfile},	#.max.vcf
		"rfile=s"	=> \$opt{rfile},	#chrM.fa
		"help"          => \$opt{help}
 	);
	die "ERROR: $! " if (!$result);
	if($opt{help})   { print $HELP; exit 0 }
	die "ERROR: $! " if (!defined($opt{mfile}));
	die "ERROR: $! " if (!defined($opt{rfile}));

	#####################################################################################

	my %diff;
	my %max;
	my %fix;

	open(IN,$opt{mfile}) or die "ERROR:$!";
        while(<IN>)
        {
		if(!/^#/)
		{
			chomp;
			my @F=split /\t/;

			# 2022/05/18; skip adjacent deletions; removed sept 10
			#next if(length($F[4])-length($F[3])<0 and $diff{$F[1]-1} and $diff{$F[1]-1}<0);

			$max{$F[1]}=$F[3];
			
			if(length($F[4])-length($F[3]))
			{
				$diff{$F[1]}=length($F[4])-length($F[3]);
			}
		}
        }
 	close(IN);

	#human
	#if($opt{rlen}==16569)
	#{
	#	if($opt{ref} eq "rCRS" or $opt{ref} eq "chrM") 
	#	{
	#		$diff{3107}=-1
	#	}
	#	elsif($opt{ref} eq "RSRS") 
	#	{
	#		$diff{523}=-1;
	#		$diff{524}=-1;
	#		$diff{3107}=-1;
	#	}
	#}
        my $MT="";
        open(IN,$opt{rfile}) or die "ERROR:$!";
        while(<IN>)
        {
                if(/^>/){}
                else
                {
                        chomp;
                        $MT.=$_;
                }
        }
	my @MT=split //,$MT;
	foreach my $i (0..@MT-1)
	{	
		$diff{$i}=-1 if(uc($MT[$i]) eq "N")
	}

	####################################################################################
	my %diff2;
	my $diff=0;
	my @pos=(sort {$a<=>$b} keys %diff);
	foreach my $pos (@pos)
	{
		$diff+=$diff{$pos};
		#$diff2{$pos+$diff+1}=$diff;  # +1: Sept 11
		$diff2{$pos+1}=$diff;	      # 2024/01/20
	}
	my @pos2=(sort {$a<=>$b} keys %diff2);

	###################################################

	while(<>)
	{
		if(/^$/)
		{
			next;
		}
		elsif(/^##reference/ or /^##contig/)
		{
			print "##reference=file://$opt{rfile}\n";
			print "##contig=<ID=$opt{ref},length=$opt{rlen}>\n";
		}
		elsif(/^#/)
		{
			print;
			next;
		}
		else
		{
			chomp;
			my @F=split /\t/;
			my $diff2=0;
			my $i=0;
			while($i<@pos2 and $F[1]>$pos2[$i])
			{
				$diff2=$diff2{$pos2[$i]};
				$i++
			}

			$F[0]=$opt{ref};
			$F[1]-=$diff2;

			my $F3=substr($MT,$F[1]-1,length($F[3]));
			
			if($F3 ne $F[3]) 
			{  
				if($F3 eq $F[4])
				{
					@F[3,4]=($F3,$F[3]); 
				}
				else
				{
					@F[3]=($F3);
				}
			}
			
		
			print join "\t",@F;
                        print "\n";
		}
	}
}

