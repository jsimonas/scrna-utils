#!/usr/bin/env python3
import os
import pysam
import argparse
import pandas as pd

# script to tag BAM generated by SEQC

def add_tags(bam, outbam, threads):
    
    sbam = bam+'sorted'
    pysam.sort("-o", sbam, bam)
    
    pysam.index(sbam)
    
    inp = pysam.AlignmentFile(sbam, 'rb', threads = threads)
    out = pysam.AlignmentFile(outbam, 'wb', template = inp, threads = threads)
        
    for read in inp.fetch():
        qname = read.query_name            
        cb = qname.split(':')[1]
        ub = qname.split(':')[2]
        read.set_tag(tag = 'CB', value = cb, value_type = 'Z')
        read.set_tag(tag = 'UB', value = ub, value_type = 'Z')
        out.write(read)
    inp.close()
    out.close()
    
    os.remove(sbam)
    os.remove(sbam+'.bai')

def main():
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument('--bam', type=str, metavar='FILENAME',
                        help='path to input BAM file')
    parser.add_argument('--out', type=str, metavar='FILENAME',
                        help='path to output bam file')
    parser.add_argument('--n', type=int, default = 4,
                        help='number of threads to use')
    
    args = parser.parse_args()

    add_tags(bam = args.bam, outbam = args.out, threads = args.n)

if __name__ == "__main__":
    main()
