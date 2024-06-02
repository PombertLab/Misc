#!/usr/bin/env python3

name = "FastqUpdate.py"
version = '0.2a'
updated = '2024-06-02'

from Bio import SeqIO
import argparse
from os import path 
from sys import exit,argv

#########################################################################
### Command line options
#########################################################################

usage = f"""
NAME        {name}
VERSION     {version}
UPDATED     {updated}
SYNOPSIS    Convert sequence file formats using biopython

REQS        python3, biopython

COMMAND     {name} \\
              --input raw_R1.fastq \\
              --encode fastq-illumina

-i (--input)    Input sequence file
-e (--encode)   Input sequence encoding
-n (--nencode)  Output sequence encoding [default = fastq]
-o (--output)   Output file name [default = <InputFileName>_converted.<nencode>]
-f (--formats)  Show supported formats
-v (--version)  Show script version
"""

formats = f"""
Supported Sequence Encoding Types
---------------------------------
r = read, w = write

r-  abi
r-  abit-trim
r-  ace
r-  cif-atom
r-  cif-seqres
rw  clustal
rw  embl
rw  fasta
rw  fasta-2line
rw  fastq-sanger/fastq
rw  fastq-solexa
rw  fastq-illumina
r-  gck
rw  genebank/gb
r-  ig
rw  imgt
rw  nexus
r-  pbd-seqres
r-  pdv-aom
rw  phd
rw  phylip
rw  pir
rw  seqxml
rw  sff
r-  sff-trim
r-  snapgene
rw  stockholm
r-  swiss
rw  tab
rw  qual 
r-  uniprot-xml 
rw  xda
"""

# Print custom message if argv is empty
if len(argv) <= 1:
    print(usage)
    exit(0)

## Setting up default variable
nencode = "fastq"

## Set up command line parser + command line flags
parser = argparse.ArgumentParser(usage=usage)
args = parser.add_argument("-i","--input")
args = parser.add_argument("-e","--encode")
args = parser.add_argument("-o","--output")
args = parser.add_argument("-n","--nencode")
args = parser.add_argument("-f","--formats", action='store_true')
args = parser.add_argument("-v","--version", action='store_true')
args = parser.parse_args()

## Parse command line objects by <parser_object>.<flag_word>
in_file = args.input
encode = args.encode
out_file = args.output
nencode = args.nencode
sformats = args.formats
scversion = args.version

#########################################################################
### Version
#########################################################################

if scversion:
    print ("")
    print (f"Script:     {name}")
    print (f"Version:    {version}")
    print (f"Updated:    {updated}\n")
    exit(0)

#########################################################################
### Supported formats
#########################################################################

if sformats:
    print (formats)
    exit(0)

#########################################################################
### File conversion
#########################################################################

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
