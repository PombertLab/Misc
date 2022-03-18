#!/usr/bin/perl
## Pombert Lab 2022
my $name = 'gc_parser.pl';
my $version = 0.1;
my $updated = '2022-03-18';

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path);

my $usage = <<"USAGE";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Parses contigs by their GC content

COMMAND		${name} \\
		  -f *.fasta \\
		  -o GC_parsed \\
		  -u 55 \\
		  -l 30

OPTIONS:
-f (--fasta)	FASTA file(s) to parse
-o (--outdir)	Output directory [Default: ./]
-u (--upper)	Upper GC (%) cutoff value [Default: 50]
-l (--lower)	Lower GC (%) cutoff value [Default: 35]
USAGE
die "\n$usage\n" unless @ARGV;

my @fasta;
my $outdir = './';
my $upper = 50;
my $lower = 35;
GetOptions(
	'f|fasta=s@{1,}' => \@fasta,
	'o|outdir=s' => \$outdir,
	'u|upper=i' => \$upper,
	'l|lower=i' => \$lower
);

## Checking for output directory
unless (-d $outdir) {
	make_path( $outdir, { mode => 0755 } )  or die "Can't create $outdir: $!\n";
}
open METRICS, ">", "$outdir/gc_metrics.tsv" or die "Can't create $outdir/gc_metrics.tsv: $!\n";
print METRICS "# Fasta file\tContig\tLength (nt)\tAT (%)\tGC (%)\n";

## Working on fasta file(s)
while (my $fasta = shift @fasta){

	my %sequences;
	my $locus;

	open FASTA, "<", $fasta or die "Can't open $fasta: $!\n";
	my ($basename) = fileparse($fasta);
	$basename =~ s/\.\w+//;
	open UPPER, ">", "$outdir/$basename.upper.fasta" or die "Can't create $basename.upper.fasta: $!\n";
	open LOWER, ">", "$outdir/$basename.lower.fasta" or die "Can't create $basename.lower.fasta: $!\n";
	open MIDR, ">", "$outdir/$basename.midrange.fasta" or die "Can't create $basename.midrange.fasta: $!\n";

	while (my $line = <FASTA>){
		chomp $line;
		if ($line =~ /^>(.*)$/){
			$locus = $1;
		}
		else {
			$sequences{$locus} .= $line;
		}
	}

	foreach my $contig (sort (keys %sequences)){

		my $sequence = $sequences{$contig};
		my $seq_length = length ($sequence);

		## Caculating nucleotide composition
		my $a_count = $sequence =~ tr/Aa//;
		my $t_count = $sequence =~ tr/Tt//;
		my $g_count = $sequence =~ tr/Gg//;
		my $c_count = $sequence =~ tr/Cc//;

		my $gc_percentage = (($g_count + $c_count)/$seq_length)*100;
		$gc_percentage = sprintf("%.2f", $gc_percentage);

		my $at_percentage = (($a_count + $t_count)/$seq_length)*100;
		$at_percentage = sprintf("%.2f", $at_percentage);

		print METRICS $fasta."\t".$contig."\t".$seq_length."\t".$at_percentage."\t".$gc_percentage."\n";

		## Writing sequences to ouput files
		if ($gc_percentage >= $upper){
			sequence(\*UPPER, $contig, $sequence);
		}
		elsif ($gc_percentage < $lower){
			sequence(\*LOWER, $contig, $sequence);
		}
		else {
			sequence(\*MIDR, $contig, $sequence);
		}

	}
}

### Subroutines
sub sequence {

	my $fh = $_[0];
	my $header = $_[1];
	my $seq = $_[2];

	print $fh ">$header\n";
	my @SEQUENCE = unpack ("(A60)*", $seq);
	while (my $seq = shift@SEQUENCE){
		print $fh "$seq\n";
	}

}

