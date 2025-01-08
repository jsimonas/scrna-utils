#!/usr/bin/env python3
import os
import pysam
import argparse
import pandas as pd
import numpy as np


def count_mapped_features(bam, threads):
    """
    Count occurrences of sF tag events translated into meaningful categories per cell barcode (CB),
    while taking into account unique molecular barcodes (UB) and gene names (GN).

    Args:
        bam (str): Input BAM file path.
        threads (int): Number of threads to use for reading the BAM file.

    Returns:
        pd.DataFrame: A DataFrame with columns for CB and all counted sF events.
    """

    # Index BAM file (if not already indexed)
    pysam.index(bam)

    # Open the BAM file
    inp = pysam.AlignmentFile(bam, 'rb', threads=threads)

    # Initialize dictionary to store unique UB and GN combinations per CB
    ub_set = {}

    # Initialize dictionary to store counts of translated sF events per CB
    sf_counts = {}

    # Translation of sF tags to meaningful categories
    translate_tags = {
        'B:i,1,1': 'fully_exonic_sense',
        'B:i,2,0': 'fully_exonic_antisense',
        'B:i,3,1': 'mainly_exonic_sense',
        'B:i,4,0': 'mainly_exonic_antisense',
        'B:i,5,1': 'mainly_intronic_sense',
        'B:i,6,0': 'mainly_intronic_antisense',
        'B:i,7,0': 'mainly_intergenic_nonsense',
        'NaN': 'unmapped'
    }

    # Parse BAM file
    for read in inp.fetch():
        # Use only primary reads
        if not read.is_supplementary and not read.is_secondary:
            # Get tags
            cb = read.get_tag('CB')
            ub = read.get_tag('UB')
            gn = read.get_tag('GN')
            sf = read.get_tag('sF') if read.has_tag('sF') else 'NaN'

            # Translate the sF tag
            sf_category = translate_tags.get(sf, 'unknown')

            # Initialize CB entry if not present
            if cb not in ub_set:
                ub_set[cb] = set()
                sf_counts[cb] = {key: 0 for key in translate_tags.values()}

            # Add UB and GN combination to the set
            if (ub, gn) not in ub_set[cb]:
                ub_set[cb].add((ub, gn))

                # Increment the count for the corresponding sF category
                sf_counts[cb][sf_category] += 1

    inp.close()

    # convert dict to dataframe
    data = []
    for cb, counts in sf_counts.items():
        row = {'barcode': cb, **counts}
        data.append(row)

    df = pd.DataFrame(data)
    df.fillna(0, inplace=True)  # Ensure no NaN values in counts

    df.to_csv(
        path_or_buf = out + '/feature_counts.tsv',
        index= False,
        sep = '\t'
    )


def main():
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument(
        '--bam', type=str, metavar='FILENAME',
        help='path to input BAM file'
    )
    parser.add_argument(
        '--out', type=str, metavar='PATH',
        help='path to output directory'
    )
    parser.add_argument(
        '--n', type=int, default = 4,
        help='number of threads to use'
    )
    
    args = parser.parse_args()

    downsample_by_cov(
        bam = args.bam,
        out = args.out,
        threads = args.n
    )

if __name__ == "__main__":
    main()
