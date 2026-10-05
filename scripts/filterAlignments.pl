#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $MINID=0;   # 0.95
my $MINPC=0;   # 0.95;
my $MINLEN=0;  # 6000;

GetOptions(
    'minid=f'  => \$MINID,
    'minpc=f'  => \$MINPC,
    'minlen=i' => \$MINLEN,
) or die "Usage: $0 [--minid F] [--minpc F] [--minlen N]\n";

while (<STDIN>) {
    chomp;

    # Keep SAM header
    if (/^\@/) {
        print "$_\n";
        next;
    }

    my @F = split /\t/;

    my $cigar = $F[5];
    my $seq   = $F[9];

    next if $cigar eq '*';
    next if $seq eq '*';

    my $read_length    = length($seq);
    my $aligned_length = 0;
    my $matches        = 0;
    my $mismatches     = 0;

    while ($cigar =~ /(\d+)([MIDNSHP=X])/g) {
      my ($length, $op) = ($1, $2);

      $aligned_length += $length if $op =~ /^[=XM]$/;
      $matches        += $length if $op =~ /^[=M]$/;
      $mismatches     += $length if $op =~ /^[X]$/;
    }

    next if $aligned_length < $MINLEN;

    my $fraction = $aligned_length / $read_length;
    next if $fraction < $MINPC;

    my $identity_count = $matches + $mismatches;
    my $identity = $matches / $identity_count;
    next if $identity < $MINID;
 
    print "$_\n";
}
