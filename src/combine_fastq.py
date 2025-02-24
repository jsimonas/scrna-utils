#!/usr/bin/env python3
import argparse
import pysam
import gzip

def combine_fastq(r1_file, r2_file, out_file):
    
    # for each read in the R1 and R2 files, extract the first 8 bases
    # (and their quality scores) from the R1 record and append them to the R2.
    with pysam.FastxFile(r1_file) as f1, pysam.FastxFile(r2_file) as f2, gzip.open(out_file, 'wt') as out:
        for r1, r2 in zip(f1, f2):
            # extract the first 8 bases (and qualities) from the R1
            r1_prefix = r1.sequence[:8]
            r1_prefix_qual = r1.quality[:8] if r1.quality is not None else ''
            
            # append the R1 prefix to the R2 record
            combined_seq = r2.sequence + r1_prefix
            combined_qual = (r2.quality if r2.quality is not None else '') + r1_prefix_qual
            
            # build the header using R2 information
            header = "@" + r2.name
            if r2.comment:
                header += " " + r2.comment
            
            out.write(f"{header}\n{combined_seq}\n+\n{combined_qual}\n")

def main():
    parser = argparse.ArgumentParser(
        description="Merge a pair of FASTQ records by appending the first 8 bases from R1 to R2."
    )
    parser.add_argument('--r1', type=str, required=True,
                        help='Path to the R1 FASTQ file (gzipped)')
    parser.add_argument('--r2', type=str, required=True,
                        help='Path to the R2 FASTQ file (gzipped)')
    parser.add_argument('--out', type=str, required=True,
                        help='Path to the output FASTQ file (gzipped)')
    parser.add_argument('--threads', type=int, default=4,
                        help='Number of threads to use (not implemented)')
    
    args = parser.parse_args()
    
    combine_fastq(args.r1, args.r2, args.out)

if __name__ == "__main__":
    main()
