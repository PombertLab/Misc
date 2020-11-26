#!/usr/bin/perl
## Pombert JF, Illinois Tech - 2019
my $version = '0.1';
my $name = 'keep_longuest_reads.pl';

use strict; use warnings; use Getopt::Long qw(GetOptions);

## Usage definition
my $usage = <<"USAGE";

NAME		$name
VERSION		$version
SYNOPSIS	Parses FASTQ file to keep the longuest reads.
		Useful for large Nanopore or PacBio FASTQ datasets.
EXAMPLE1		keep_longuest_reads.pl -i file.fastq -o parsed.10k.fastq -m 10000
EXAMPLE2		keep_longuest_reads.pl -i file.fastq -o parsed.100X.fastq -d 100 -s 3000000 

OPTIONS
-i (--input)	Input file in FASTQ format
-o (--output)	Output file name
-m (--minimum)	Minimum read length to keep
-d (--depth)	Desired sequencing depth (requires estimated genome size: -s)
-s (--size)	Expected genome size

NOTE: -m and -d are mutually exclusive
USAGE

die "$usage\n" unless@ARGV;

my $fastq;
my $output;
my $min;
my $depth;
my $genome_size;
GetOptions(
	'i|input=s' => \$fastq,
	'o|output=s' => \$output,
	'm|minimum=i' => \$min,
	'd|depth=i' => \$depth,
	's|size=i' => \$genome_size
);

open FASTQ1, "<$fastq";
open FASTQ2, "<$fastq";
open OUT, ">$output";

## parsing by minimum length
if ($min){
	my %reads; my $count = 0; my $read;
	while (my $line = <FASTQ1>){
		chomp $line;
		if ($count == 0){$read = $line; $reads{$read}[0] = $line; $count++;} ## Read name
		elsif ($count == 1){$reads{$read}[1] = $line; $count++;} ## Read sequence
		elsif ($count == 2){$count++;} ## Skipping + sign
		elsif ($count == 3){
			$reads{$read}[2] = $line; ## Quality score
			if (length($reads{$read}[1]) >= $min){
				print OUT "$reads{$read}[0]\n";
				print OUT "$reads{$read}[1]\n";
				print OUT '+'."\n";
				print OUT "$reads{$read}[2]\n";
			}
			$count = 0; %reads = (); ## Clearing read db to minimize memory usage
		}
	}
}
elsif ($depth){
	my @lengths; my $count = 0;
	## Pass 1 - Calculating metrics; doing two independent passes to avoid loading the full FASTQ file in memory
	while (my $line = <FASTQ1>){
		chomp $line;
		if ($count == 0){$count++;}
		elsif ($count == 1){
			my $len = sprintf("%09d", length($line)); ## Adding leading zeroes to help sort array
			push (@lengths, $len);
			$count++;
		}
		elsif($count == 2){$count++;}
		elsif($count == 3){$count=0;}
	}
	@lengths = sort @lengths; ## sort by size
	@lengths = reverse @lengths; ## from largest to smallest
	my $len_threshold; my $sum;
	while (my $x = shift @lengths){
		if ($sum <= ($depth*$genome_size)){$len_threshold = $x;}
		$sum += $x;
	}
	$len_threshold = sprintf("%0d", $len_threshold);
	print "\nMinimum read size for ${depth}X sequencing depth at estimated genome size of $genome_size bp = $len_threshold bp\n";
	print "Saving reads of at least $len_threshold bp to $output\n\n";
	## Pass 2 - Parsing reads
	$count = 0; my $read; my %reads; 
	while (my $line = <FASTQ2>){
		chomp $line;
		if ($count == 0){$read = $line; $reads{$read}[0] = $line; $count++;} ## Read name
		elsif ($count == 1){$reads{$read}[1] = $line; $count++;} ## Read sequence
		elsif ($count == 2){$count++;} ## Skipping + sign
		elsif ($count == 3){
			$reads{$read}[2] = $line; ## Quality score
			if (length($reads{$read}[1]) >= $len_threshold){
				print OUT "$reads{$read}[0]\n";
				print OUT "$reads{$read}[1]\n";
				print OUT '+'."\n";
				print OUT "$reads{$read}[2]\n";
			}
			$count = 0; %reads = (); ## Clearing read db to minimize memory usage
		}
	}
}
