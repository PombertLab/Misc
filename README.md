## Misc
### Synopses

`FastqUpdate.py`: 
Updates old Illumina FASTQ files to the new Illumina/Sanger 1.9 FASTQ encoding format. Requires biopython.

`gc_parser.pl`:
Parses multifasta files by their GC content. Useful to parse out contaminants from genome assemblies.

`keep_longest_reads.pl`: Calculates metrics for FASTQ file(s) and/or parses them to keep the longest reads either by minimum size or by desired sequencing depth (useful for large Nanopore or PacBio datasets).

`read_len_plot.py`: Plots the read length distribution for a given FASTQ dataset with matplotlib.

`run_pilon.pl`: Runs [Pilon](http://software.broadinstitute.org/software/pilon/) read correction in paired-end mode for X iterations (stops automatically if Pilon no longer makes changes to the consensus)

`runTaxonomizedBLAST.pl`: Runs taxonomized [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) searches, and returns the outfmt 6 format with columns staxids, sscinames, sskingdoms, and sblastnames.

`parseTaxonomizedBLAST.pl`: Parses the content of [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) searches performed with `runTaxonomizedBLAST.pl`.