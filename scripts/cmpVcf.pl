#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that compares a REF to a QRY vcf file and generates a confusion matrix

        EXAMPLE:
                $0 I.vcf J.vcf 

~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
        my (%opt,%h0,%h1);

        # validate input parameters
        my $result = GetOptions(
		"sm"     =>     \$opt{sm},
		"help"	 =>	\$opt{help},
		"snp"	 => 	\$opt{snp},
		"ins"	 => 	\$opt{ins},
                "del"    => 	\$opt{del},
		"noh"	 =>	\$opt{noh},
		"het"	 =>	\$opt{het},
                "hom"    =>     \$opt{hom},
		"cds"	 =>	\$opt{cds},
		"gene"   =>	\$opt{gene},

		"pos"	=>	\$opt{pos}
	);
	$opt{pos}=1;
	#$opt{noh}=1;

	if(!$result)            { die "ERROR: $! "}
        if($opt{help})          { print $HELP; exit 0 }
        if(@ARGV<2)             { die "ERROR: insufficient number of arguments"}

        #########################################

        open(IN,$ARGV[1]) or die("ERROR: Cannot open input file".$!) ;
        while(<IN>)
        {
                next if(/^$/ or /^#/);

		chomp;
                my @F=split /\t/;
                die "ERROR $_" if(@F<5);
		next if($F[3]=~/N/);

                if($opt{snp})   { next if(length($F[3]) ne length($F[4])) }
                if($opt{ins} )  { next if(length($F[3]) ge length($F[4])) }
                if($opt{del})   { next if(length($F[3]) le length($F[4])) }
		next if($opt{cds}  and $_!~/CDS=/);
		next if($opt{gene} and $_!~(/CDS=/ or /TRN=/ or /RNR=/));
		#next if($opt{noh}  and /Homopolymer/i);
		next if($opt{noh}  and (296<$F[1] and $F[1]<318 or 493<$F[1] and $F[1]<502 or 511<$F[1] and $F[1]<524 or 538<$F[1] and $F[1]<545 or 561<$F[1] and $F[1]<574 or 954<$F[1] and $F[1]<965 or 5888<$F[1] and $F[1]<5895 or 8269<$F[1] and $F[1]<8288 or 16178<$F[1] and $F[1]<16193));

		next if($opt{het} and !(/\tAF=0/ or /;AF=0/ or /:0\.\d+$/));
                next if($opt{hom} and (/\tAF=0/ or /;AF=0/ or /:0\.\d+$/));


		my $SM="";

		if($opt{sm})
		{
			if(@F>8 and $F[8] eq "SM")                                   { $SM=$F[9] }
			elsif(@F>7 and ($F[7]=~/SM=(\S+?);/ or $F[7]=~/SM=(\S+?)$/)) { $SM=$1 }
			else						             { die "ERROR: $_" }
		}

		my $AF=1;
		$AF=$1 if(/AF=(0\.\d+)/ or /.+:(1)$/ or /.+:(0\.\d+)/);
		($F[3],$F[4])=(uc($F[3]),uc($F[4]));

		if($opt{pos}) { $h1{"$F[0] $F[1] $SM"}=$AF; }
                else          { $h1{"$F[0] $F[1] $F[3] $F[4] $SM"}=$AF; }
        }
	close(IN);
        #last unless(%h1);

        #########################################

        open(IN,$ARGV[0]) or die("ERROR: Cannot open input file".$!) ;
        while(<IN>)
        {
                next if(/^$/ or /^#/);
                chomp;
                my @F=split /\t/;
                die "ERROR $_" if(@F<5);
		next if($F[3]=~/N/);

                if($opt{snp})   { next if(length($F[3]) ne length($F[4])) }
                if($opt{ins})   { next if(length($F[3]) ge length($F[4])) }
                if($opt{del})   { next if(length($F[3]) le length($F[4])) }
                next if($opt{cds}  and $_!~/CDS=/);
                next if($opt{gene} and $_!~(/CDS=/ or /TRN=/ or /RNR=/));
		#next if($opt{noh}  and /Homopolymer/i);
		next if($opt{noh}  and (296<$F[1] and $F[1]<318 or 493<$F[1] and $F[1]<502 or 511<$F[1] and $F[1]<524 or 538<$F[1] and $F[1]<545 or 561<$F[1] and $F[1]<574 or 954<$F[1] and $F[1]<965 or 5888<$F[1] and $F[1]<5895 or 8269<$F[1] and $F[1]<8288 or 16178<$F[1] and $F[1]<16193));

                next if($opt{het} and !(/\tAF=0/ or /;AF=0/ or /:0\.\d+$/));
                next if($opt{hom} and (/\tAF=0/ or /;AF=0/ or /:0\.\d+$/));
		
                my $SM="";

                if($opt{sm})
                {
                        if(@F>8 and $F[8] eq "SM")                                   { $SM=$F[9] }
                        elsif(@F>7 and ($F[7]=~/SM=(\S+?);/ or $F[7]=~/SM=(\S+?)$/)) { $SM=$1 }
                        else                                                         { die "ERROR: $_" }
                }

                my $AF=1;
                $AF=$1 if(/AF=(0\.\d+)/ or /.+:(1)$/ or /.+:(0\.\d+)/);
                ($F[3],$F[4])=(uc($F[3]),uc($F[4]));

                if($opt{pos}) {	$h0{"$F[0] $F[1] $SM"}=$AF; }
                else          { $h0{"$F[0] $F[1] $F[3] $F[4] $SM"}=$AF; }
        }
        close(IN);

	#########################################
        my ($TP,$FP,$FN,$S,$P,$F1)=(0,0,0,".",".",".");

	foreach my $k (keys %h0)
	{
		if($h1{$k}) 
		{ 
			$TP++
		}
		else	    
		{ 
			$FN++
		}
	}

	foreach my $k (keys %h1)
        {
	        if(!$h0{$k}) 
		{ 
			$FP++
        	}
	}

	$S=$TP/($TP+$FN)    if($TP+$FN);
	$P=$TP/($TP+$FP)    if($TP+$FP);
	$F1=2*$S*$P/($S+$P) if($S=~/\d/ and $P=~/\d/ and $S+$P);

	$S=int($S*100+.5)/100   unless($S eq ".");
	$P=int($P*100+.5)/100   unless($P eq ".");
	$F1=int($F1*100+.5)/100 unless($F1 eq ".");

	my $snv=".";
	if($opt{snp})      { $snv="snp"}
	elsif($opt{ins}) { $snv="ins"}        
        elsif($opt{del}) { $snv="del"}	
	
	if($opt{noh}) { $opt{noh}=1}
	else {$opt{noh}="."}

        print join "\t",("qry","snv","TP","FN","FP","S","P","F1"); print "\n";
        print join "\t",($ARGV[1],$snv,$TP,$FN,$FP,$S,$P,$F1); print "\n";

	exit 0;
}
