#!/usr/bin/env python

# Description: Convert build coordinates from hg19 to hg38 with CrossMap.py (pip install crossmap) https://crossmap.readthedocs.io/
# Usage: python convert_build_coordinates.py -i input_file -o output_file
# Input file: tab-delimited file with at least two columns: CHROM and GENPOS
# Output file: tab-delimited file with two additional columns: chr_hg38 and pos_hg38

# Note: The output file will contain duplicated coordinates if the input file contains duplicated coordinates
# Note: The output file will remove unmapped coordinates if the input file contains unmapped coordinates
# Note: A liftover conversion file needs to be specified with --conversion; Download from https://crossmap.readthedocs.io/en/latest/#chain-file

# By: Chang Lu

import pandas as pd
import subprocess
import argparse

def main():
    # Parse command line arguments 
    parser = argparse.ArgumentParser(description='Convert build coordinates')
    parser.add_argument('-i', '--input', required=True, help='Input file')
    parser.add_argument('-o', '--output', required=True, help='Output file')
    parser.add_argument('--chr', default='CHROM', help='Chromosome prefix')
    parser.add_argument('--pos', default='GENPOS', help='Position column')
    parser.add_argument('--conversion', default='/Users/clu/work/genomic_resource/hg19ToHg38.over.chain.gz', help='Build conversion file')
    args = parser.parse_args()

    # Read input file
    df = pd.read_csv(args.input, sep='\t')

    # Position + 1
    pos1 = df[args.pos] + 1
    
    # Prepare BED file
    beddf = df[[args.chr, args.pos]].assign(pos1=pos1)
    # remove duplicated coordinates
    beddf = beddf.drop_duplicates()

    # Write coordinates into a temporary file
    beddf.to_csv('temp.bed', sep='\t', index=False, header=False)
    
    # Convert build coordinates
    cmd = 'CrossMap.py bed {} temp.bed > temp2.bed'.format(args.conversion)
    subprocess.call(cmd, shell=True)

    # Check unmap coordinates
    unmap_outpath = '{}_unmap.bed'.format(args.output)
    cmd = 'grep Unmap temp2.bed > {}'.format(unmap_outpath)
    subprocess.call(cmd, shell=True)
    # Check if temp_unmap.bed is empty with python
    if open(unmap_outpath).read() == '':
        df_unmap = None
        print('All coordinates are mapped')
    else:
        # Unmapped coordinates
        df_unmap = pd.read_csv(unmap_outpath, sep='\t', header=None)
        df_unmap = df_unmap.assign(hg19=['{}:{}'.format(x, y) for x, y in zip(df_unmap[0], df_unmap[1])])
        print('Unmapped coordinates: {}'.format(df_unmap.shape[0]))

    # Read converted coordinates
    cmd = 'grep -v Unmap temp2.bed > temp_map.bed'
    subprocess.call(cmd, shell=True)
    df_map = pd.read_csv('temp_map.bed', sep='\t', header=None)

    # Merge Mapped coordinates
    df_map = df_map.assign(hg19=['{}:{}'.format(x, y) for x, y in zip(df_map[0], df_map[1])])
    df_map = df_map.set_index('hg19')
    df_map = df_map.assign(chr_hg38=df_map[4], pos_hg38=df_map[5])
    df_map = df_map[['chr_hg38', 'pos_hg38']]
    df_map = df_map.astype({'pos_hg38': 'int32'})

    # All coordinates
    df = df.assign(hg19=['{}:{}'.format(x, y) for x, y in zip(df[args.chr], df[args.pos])])
    df = df.set_index('hg19')

    # Remove unmapped coordinates
    if df_unmap is not None:
        df_unmap = df_unmap.set_index('hg19')
        df = df.drop(df_unmap.index)

    # Join by index and drop left
    df = df.join(df_map, how='left')
    
    # Write output file
    df.to_csv(args.output, sep='\t', index=False)

    # Remove temporary files
    cmd = 'rm temp.bed temp2.bed temp_map.bed'
    subprocess.call(cmd, shell=True)

if __name__ == '__main__':
    main()