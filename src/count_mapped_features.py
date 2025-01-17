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

    # sF tag categories based on the first digit
    translate_tags = {
        '1': 'fully_exonic_se',
        '2': 'fully_exonic_as',
        '3': 'mainly_exonic_se',
        '4': 'mainly_exonic_as',
        '5': 'mainly_intronic_se',
        '6': 'mainly_intronic_as',
        '7': 'unique_intergenic_ns',
        '-1': 'multi_intergenic_ns',
        'unknown': 'unknown'
    }

    # parse BAM
    for read in inp.fetch(until_eof=True):
        # use only primary reads
        if not read.is_supplementary and not read.is_secondary:
            # get tag
            cb = read.get_tag('CB')

            # init CB entry if not present
            if cb not in sf_counts:
                sf_counts[cb] = {category: 0 for category in translate_tags.values()}
                sf_counts[cb]['unknown'] = 0
                sf_counts[cb]['unmapped_reads'] = 0
                total_reads[cb] = 0

            # check if the read is unmapped
            if read.is_unmapped:
                sf_counts[cb]['unmapped_reads'] += 1
                total_reads[cb] += 1
                continue
            
            # process mapped reads
            sf = ','.join(map(str, read.get_tag('sF'))) if read.has_tag('sF') else 'NaN'
            sf_key = sf.split(',')[0] if sf != 'NaN' else 'NaN'
            sf_category = translate_tags.get(sf_key, 'unknown')
    
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
