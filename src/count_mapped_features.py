#!/usr/bin/env python3
import os
import pysam
import argparse
import pandas as pd
import numpy as np

def count_mapped_features(bam, out, threads):
    """
    Count occurrences of sF tag events translated into meaningful categories per cell barcode (CB),
    and add a total_reads column for each CB.

    Args:
        bam (str): Input BAM file path.
        threads (int): Number of threads to use for reading the BAM file.

    Returns:
        Tab separated DataFrame with columns for cell barcodes and all counted sF events.
    """

    # index BAM
    pysam.index(bam)

    # open BAM
    inp = pysam.AlignmentFile(bam, 'rb', threads=threads)

    # init dicts
    sf_counts = {}
    total_reads = {}

    # sF tag categories
    translate_tags = {
        '1,1': 'fully_exonic_se',
        '2,0': 'fully_exonic_as',
        '3,1': 'mainly_exonic_se',
        '4,0': 'mainly_exonic_as',
        '5,1': 'mainly_intronic_se',
        '6,0': 'mainly_intronic_as',
        '7,0': 'unique_intergenic_ns',
        '-1,-1': 'multi_intergenic_ns',
        'NaN': 'unmapped',
        'unknown': 'unknown'
    }

    # parse BAM
    for read in inp.fetch():
        # use only primary reads
        if not read.is_supplementary and not read.is_secondary:
            # get tags
            cb = read.get_tag('CB')
            sf = ','.join(map(str, read.get_tag('sF'))) if read.has_tag('sF') else 'NaN'
            
            # translate the sF tag
            sf_category = translate_tags.get(sf, 'unknown')

            # init CB entry if not present
            if cb not in sf_counts:
                sf_counts[cb] = {key: 0 for key in translate_tags.values()}
                total_reads[cb] = 0

            # increment the count for the corresponding sF category
            sf_counts[cb][sf_category] += 1

            # increment total_reads for the CB
            total_reads[cb] += 1

    inp.close()

    # convert dict to dataframe
    data = []
    for cb, counts in sf_counts.items():
        row = {'barcode': cb, 'total_reads': total_reads[cb], **counts}
        data.append(row)

    df = pd.DataFrame(data)
    df.fillna(0, inplace=True)

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

    count_mapped_features(
        bam = args.bam,
        out = args.out,
        threads = args.n
    )

if __name__ == "__main__":
    main()
