#!/usr/bin/python
from Bio import SeqIO
import argparse
from os import path 
from sys import exit,argv

name = "FastqUpdate.py"
version = '0.2'
updated = '2021-03-21'

usage = f"""
NAME		{name}
VERSION		{version}
UPDATED		{updated}
SYNOPSIS	The purpse of this script is to convert sequence files of one encoding
			into sequence files of another encoding type

COMMAND		{name} \\
			--input raw_R1.fastq \\
			--encode fastq-illumina

-i | --input	Input sequence file
-e | --encode	Input sequence encoding
-n | --nencode	Output sequence encoding [default = fastq]
-o | --output	Output file name [default = <InputFileName>_converted.<nencode>]

Supported Sequence Encoding Types
---------------------------------
KEY: r|ead w|rite

r-|abi
r-|abit-trim
r-|ace
r-|cif-atom
r-|cif-seqres
rw|clustal
rw|embl
rw|fasta
rw|fasta-2line
rw|fastq-sanger/fastq
rw|fastq-solexa
rw|fastq-illumina
r-|gck
rw|genebank/gb
r-|ig
rw|imgt
rw|nexus
r-|pbd-seqres
r-|pdv-aom
rw|phd
rw|phylip
rw|pir
rw|seqxml
rw|sff
r-|sff-trim
r-|snapgene
rw|stockholm
r-|swiss
rw|tab
rw|qual 
r-|uniprot-xml 
rw|xda

"""

# Die usage unless argv is populated
if len(argv) == 1:
    print(usage)
    exit()

## Setting up default variable
nencode = "fastq"

## Set up command line parser
parser = argparse.ArgumentParser(usage=usage)

## Create command line flag
args = parser.add_argument("-i","--input",required=True)
args = parser.add_argument("-e","--encode",required=True)
args = parser.add_argument("-o","--output")
args = parser.add_argument("-n","--nencode")

## Create parser object
args = parser.parse_args()

## Parse command line objects by <parser_object>.<flag_word>
in_file = args.input
encode = args.encode
out_file = args.output
nencode = args.nencode

## Parsing name and extension of input file provided
base_name, file_name = path.split(in_file)
file_name, file_ext = path.splitext(file_name)

## If output file not given, create default from input file
if not out_file:
	out_file = file_name + "_converted" + nencode

## If directory of file is current working folder, make note
if not base_name:
    base_name = "."

SeqIO.convert(in_file,encode,out_file,nencode)
