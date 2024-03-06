#!/usr/bin/perl -w
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that generates SNV count statistics
	H:Homoplasmies	
	h:heteroplasmies
	S:Homoplasmies(no INDELs)
	S:heteroplasmies(no INDELs)
	I:Homoplasmies(INDELs)
        i:heteroplasmies(INDELs)

	*p:same as aboble, non-homopolymeric regions

	A:all

        EXAMPLE:
                cat I.tab | $0 > O.summary

		Run	 H	h	S	s	I	i	Hp	hp	Sp	sp	Ip	ip	A
		HG00438	 35	2	32	0	3	2	31	0	30	0	1	0	37
		HG00621	 26	4	25	3	1	1	24	2	24	2	0	0	30
		...
		
		=>

		id	count	nonZero	min	max	median	mean	sum
		H	40	40	16	84	43	46.4	1856
		h	40	38	0	6	2	2.33	93
		...

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
	my @keys=();
	my %vals;

	# validate input parameters
        my $result = GetOptions(
		"p=s"	=> \$options{prefix},
		"help"  => \$options{help}
	);
	if(!$result)             { die "ERROR: $! "}
        if($options{help})       { print $HELP; exit 0 }

	########################################

	while(<>)
	{
		chomp;
		if($.==1)
		{
			@keys=split /\t/;
			shift @keys;
		}	
		else
		{
			my @F=split /\t/;
			shift @F;
	
			foreach my $i (0..@F-1)
			{
				push @{$vals{$keys[$i]}},$F[$i];
			}		
		}		 
	}
		
	#########################################	
	
	print "$options{prefix}\t" if($options{prefix});
	print join "\t", ("id","count","nonZero","min","max","median","mean","sum");	
	print "\n";
		
	foreach my $key (@keys)
	{	
		next unless($vals{$key});

		my @vals=sort {$a<=>$b} @{$vals{$key}};		
		my $count=scalar(@vals);
		
		my $sum=0;
		my $nonZero=0;
		foreach my $val (@vals) 
		{ 
			$sum+=$val ;
			$nonZero++ if($val);
		}
							
		my $min=$vals[0];
		my $max=$vals[-1];
	
		my $index=int(scalar(@vals)/2);
		my $median=$vals[$index];
                my $mean=int($sum*100/$count+.5)/100;
		
		print "$options{prefix}\t" if($options{prefix});
		print join "\t", ($key,$count,$nonZero,$min,$max,$median,$mean,$sum);
		print "\n"
	}		
	

	exit 0;
}
