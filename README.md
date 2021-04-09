## Misc
### Synopses

[FastqUpdate.py](https://github.com/PombertLab/Misc/blob/main/FastqUpdate.py) updates old Illumina FASTQ files to the new Illumina/Sanger 1.9 FASTQ encoding format. Requires biopython.

[keep_longest_reads.pl](https://github.com/PombertLab/Misc/blob/main/keep_longest_reads.pl) parses FASTQ files to keep the longest reads by size or desired sequencing depth. Useful for large Nanopore or PacBio FASTQ datasets.

[runTaxonomizedBLAST.pl](https://github.com/PombertLab/Misc/blob/main/runTaxonomizedBLAST.pl) runs taxonomized [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) searches, and returns the outfmt 6 format with columns staxids, sscinames, sskingdoms, and sblastnames.

[parseTaxonomizedBLAST.pl](https://github.com/PombertLab/Misc/blob/main/parseTaxonomizedBLAST.pl) parses the content of taxonomized [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) searches performed with [runTaxonomizedBLAST.pl](https://github.com/PombertLab/Misc/blob/main/runTaxonomizedBLAST.pl).

[run_pilon.pl](https://github.com/PombertLab/Misc/blob/main/run_pilon.pl) runs [Pilon](http://software.broadinstitute.org/software/pilon/) read correction in paired-end mode for X iterations (stops automatically if Pilon no longer makes changes to the consensus)