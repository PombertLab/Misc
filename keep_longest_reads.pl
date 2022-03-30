#!/usr/bin/perl
## Pombert JF, Illinois Tech - 2020
my $version = '0.7a';
my $name = 'keep_longest_reads.pl';
my $updated = '2022-03-30';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use PerlIO::gzip;

## Usage definition
my $usage = <<"USAGE";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Parses FASTQ file to keep the longest reads by size or desired sequencing depth.
		Useful for large Nanopore or PacBio FASTQ datasets.

EXAMPLE 1	${name} -i file.fastq -o parsed.10k.fastq -m 10000
EXAMPLE 2	${name} -i file.fastq -o parsed.100x.fastq -d 100 -s 3000000
EXAMPLE 3	${name} -i file.fastq -x

OPTIONS
-i (--input)	Input file in FASTQ format
-x (--metrics)	Calculate metrics only, do not create read subset
-o (--output)	Output file name
-m (--minimum)	Minimum read length to keep
-d (--depth)	Desired sequencing depth (requires estimated genome size: -s)
-s (--size)	Expected genome size

NOTE: -m and -d are mutually exclusive
USAGE
die "\n$usage\n" unless@ARGV;

my @commands = @ARGV;
my $fastq;
my $metrics;
my $output;
my $min;
my $depth;
my $genome_size;
GetOptions(
	'i|input=s' => \$fastq,
	'x|metrics'	=> \$metrics,
	'o|output=s' => \$output,
	'm|minimum=i' => \$min,
	'd|depth=i' => \$depth,
	's|size=i' => \$genome_size
);

## Input file
my $gzip = '';
if ($fastq =~ /.gz$/){ $gzip = ':gzip'; }
open FASTQ1, "<$gzip", "$fastq" or die "Can't open $fastq: $!\n";
open FASTQ2, "<$gzip", "$fastq" or die "Can't open $fastq: $!\n";

## Log file
my $logfile;
if ($metrics){ $logfile = 'metrics.log'; }
else  { $logfile = "$output.log"; }
open LOG, ">", "$logfile" or die "Can't create $logfile: $!\n";

## Output file
unless ($metrics){
	open OUT, ">", "$output" or die "Can't create $output: $!\n";
}

## Log file
my $stime = `date`; chomp $stime;
print LOG "COMMAND: $name @commands\n";
print LOG "Started on $stime\n";

## parsing by minimum length
if ($min or $metrics){

    my @lengths; my @subset; my %reads; my $count = 0; my $read;

	if ($metrics){ $min = 1; }

	while (my $line = <FASTQ1>){
		chomp $line;

		## Read name
		if ($count == 0){
			$read = $line;
			$reads{$read}[0] = $line;
			$count++;
		}

		## Read sequence
		elsif ($count == 1){ 
			$reads{$read}[1] = $line;
			## Adding leading zeroes to help sort array
			my $len = sprintf("%09d", length($line));
			push (@lengths, $len);
			$count++;
		}

		## Skipping + sign
		elsif ($count == 2){ $count++; }

		## Quality score 
		elsif ($count == 3){
			$reads{$read}[2] = $line;
			if (length($reads{$read}[1]) >= $min){
				my $keep = sprintf("%09d", length($reads{$read}[1]));
				push (@subset, $keep);
				unless ($metrics){
					print OUT "$reads{$read}[0]\n";
					print OUT "$reads{$read}[1]\n";
					print OUT '+'."\n";
					print OUT "$reads{$read}[2]\n";
				}
			}
			## Clearing read db to minimize memory usage
			$count = 0; %reads = (); 
		}
	}

	if ($fastq =~ /.gz$/){ binmode FASTQ1, ":gzip(none)"; }	

	n50($fastq, @lengths);
	unless ($metrics){
		n50($output, @subset);
	}
}

elsif ($depth){
	my @lengths;
	my $count = 0;
	## Doing two independent passes to avoid loading the full FASTQ file in memory
	## Pass 1 - Calculating metrics
	while (my $line = <FASTQ1>){
		chomp $line;
		if ($count == 0){ $count++; }
		elsif ($count == 1){
			## Adding leading zeroes to help sort array
			my $len = sprintf("%09d", length($line));
			push (@lengths, $len);
			$count++;
		}
		elsif ($count == 2){ $count++; }
		elsif ($count == 3){ $count=0; }
	}

	if ($fastq =~ /.gz$/){ binmode FASTQ1, ":gzip(none)"; }	

	@lengths = sort @lengths; ## sort by size
	@lengths = reverse @lengths; ## from largest to smallest
	my $len_threshold; my $sum;
	foreach (@lengths){
		if ($sum <= ($depth*$genome_size)){ $len_threshold = $_; }
		$sum += $_;
	}
	$len_threshold = sprintf("%0d", $len_threshold);
	print "\nMinimum read size for ${depth}X sequencing depth at estimated genome size of $genome_size bp = $len_threshold bp\n";
	print "Saving reads of at least $len_threshold bp to $output\n\n";

	## Pass 2 - Parsing reads
	$count = 0; my $read; my %reads; my @subset;
	while (my $line = <FASTQ2>){
		chomp $line;

		## Read name
		if ($count == 0){
			$read = $line;
			$reads{$read}[0] = $line;
			$count++;
		}

		## Read sequence
		elsif ($count == 1){ 
			$reads{$read}[1] = $line;
			$count++;
		}

		## Skipping + sign
		elsif ($count == 2){ $count++; }

		## Quality score
		elsif ($count == 3){
			$reads{$read}[2] = $line; 
			if (length($reads{$read}[1]) >= $len_threshold){
				my $keep = sprintf("%09d", length($reads{$read}[1]));
				push (@subset, $keep);
				print OUT "$reads{$read}[0]\n";
				print OUT "$reads{$read}[1]\n";
				print OUT '+'."\n";
				print OUT "$reads{$read}[2]\n";
			}
			## Clearing read db to minimize memory usage
			$count = 0; %reads = ();
		}
	}

	if ($fastq =~ /.gz$/){ binmode FASTQ2, ":gzip(none)"; }	

	n50($fastq, @lengths);
	n50($output, @subset);
}

my $etime = `date`; chomp $etime;
unless ($metrics){
	print LOG "Ended on: $etime\n";
}

### subroutines
sub n50{
	my @fh = (*LOG, *STDOUT);
	my $file = shift @_;
	my $num_reads = scalar @_;
	my @len = sort @_; ## sort by size
	@len = reverse @len; ## from largest to smallest

	my $nreads = commify($num_reads);
	foreach (@fh){
		print $_ "\n## Metrics for dataset $file\n\n";
		print $_ "Number of reads: $nreads\n";
	}

	## Median
	my $median;
	my $median_pos = $num_reads/2;
	if ($median_pos =~ /^\d+$/){ $median = $len[$median_pos]; }
	else {
		my $med1 = int($median_pos);
		my $med2 = $med1 + 1;
		$median = (($len[$med1] + $len[$med2])/2);
	}

	## Average
	my $sum; foreach (@len){ $sum += $_; }
	my $fsum = commify($sum);
	my $large = sprintf("%.0f", $len[0]); $large = commify($large);
	my $small = sprintf("%.0f", $len[$#len]); $small = commify($small);
	my $average = sprintf("%.0f", ($sum/$num_reads)); $average = commify($average);
	$median = sprintf("%.0f", $median); $median = commify($median);
	foreach (@fh){
		print $_ "Total number of bases: $fsum\n";
		print $_ "Longest read = $large nt\n";
		print $_ "Shortest read = $small nt\n";
		print $_ "Average read size = $average nt\n";
		print $_ "Median read size = $median nt\n";
	}

	## N50, N75, N90
	my $n50_td = $sum*0.5; my $n75_td = $sum*0.75; my $n90_td = $sum*0.9;
	my $n50; my $n75, my $n90;
	my $nsum50 = 0; my $nsum75 = 0; my $nsum90 = 0;
	foreach (@len){ $nsum50 += $_; if ($nsum50 >= $n50_td){ $n50 = $_; last; }}
	foreach (@len){ $nsum75 += $_; if ($nsum75 >= $n75_td){ $n75 = $_; last; }}
	foreach (@len){ $nsum90 += $_; if ($nsum90 >= $n75_td){ $n90 = $_; last; }}
	$n50 = sprintf ("%.0f", $n50); $n50 = commify($n50);
	$n75 = sprintf ("%.0f", $n75); $n75 = commify($n75);
	$n90 = sprintf ("%.0f", $n90); $n90 = commify($n90);
	foreach (@fh){ print $_ "N50 = $n50 nt\n"."N75 = $n75 nt\n"."N90 = $n90 nt\n"."\n"; }
}

sub commify { ## From the Perl Cookbook; O'Reilly
	my $text = reverse $_[0];
	$text =~ s/(\d{3})(?=\d)(?!\d*\.)/$1,/g;
	return scalar reverse $text;
}
