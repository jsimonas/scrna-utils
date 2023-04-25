#!/usr/bin/env python3
import os
import sys
import pysam
import argparse
import pandas as pd

# script to trim query reads in BAM

def trim_reads(bam, outbam, trim_length, threads):
    
    pysam.index(bam)
    
    inp = pysam.AlignmentFile(bam, 'rb', threads = threads)
    out = pysam.AlignmentFile(outbam, 'wb', template = inp, threads = threads)
        
    for read in inp.fetch():
        qual = read.qual[:trim_length]
        read.query_sequence = read.query_sequence[:trim_length]
        read.cigarstring = str(trim_length) + "M"
        read.qual = qual
        out.write(read)
    inp.close()
    out.close()


def main():
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument('--bam', type=str, metavar='FILENAME',
                        help='path to input BAM file')
    parser.add_argument('--out', type=str, metavar='FILENAME',
                        help='path to output bam file')
    parser.add_argument('--trim_len', type=int, default = 35,
                        help='length to trimm query read')
    parser.add_argument('--n', type=int, default = 4,
                        help='number of threads to use')
    
    args = parser.parse_args()

    trim_reads(bam = args.bam, outbam = args.out, trim_length = args.trim_len, threads = args.n)

if __name__ == "__main__":
    main()
