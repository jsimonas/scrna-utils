#!/usr/bin/env python3
import os
import sys
import pysam
import argparse
import pandas as pd

# script to tag BAM generated by SEQC

def add_tags(bam, outbam, platform, threads):
    
    sbam = bam+'sorted'
    pysam.sort("-o", sbam, bam)
    
    pysam.index(sbam)
    
    inp = pysam.AlignmentFile(sbam, 'rb', threads = threads)
    out = pysam.AlignmentFile(outbam, 'wb', template = inp, threads = threads)
        
    for read in inp.fetch():
        qname = read.query_name            
        cx = qname.split(':')[1]
        ub = qname.split(':')[2]
        if len(cx) == 16:
            cb = cx[:8]+'GAGTGATTGCTTGTGACGCCTT'+cx[8:]+ub
        elif len(cx) == 17:
            cb = cx[:9]+'GAGTGATTGCTTGTGACGCCTT'+cx[9:]+ub
        elif len(cx) == 18:
            cb = cx[:10]+'GAGTGATTGCTTGTGACGCCTT'+cx[10:]+ub
        else:
            cb = cx[:11]+'GAGTGATTGCTTGTGACGCCTT'+cx[11:]+ub
        if platform == 'indrops':
            read.set_tag(tag = 'CB', value = cb, value_type = 'Z')
        elif platform == '10x':
            read.set_tag(tag = 'CB', value = cx, value_type = 'Z')
        else:
            sys.exit('provided platform is not supported, use indrops or 10x.')
        read.set_tag(tag = 'UB', value = ub, value_type = 'Z')
        read.query_name = qname.split(';')[1]
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
    parser.add_argument('--platform', type=str, default = 'indrops',
                        help='indrops or 10x')
    parser.add_argument('--n', type=int, default = 4,
                        help='number of threads to use')
    
    args = parser.parse_args()

    add_tags(bam = args.bam, outbam = args.out, platform = args.platform, threads = args.n)

if __name__ == "__main__":
    main()
