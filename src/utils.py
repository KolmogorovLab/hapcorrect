import numpy as np
import pandas as pd
import pysam
import os

from extras import get_contigs_list
from process_bam import get_segments_coverage

def get_chromosomes_bins(bam_file, bin_size, arguments):
    bed=[]
    bam_alignment = pysam.AlignmentFile(bam_file)
    headers = bam_alignment.header
    seq_dict = headers['SQ']
    region = [''] * len(seq_dict)
    chrs = [''] * len(seq_dict)
    head, tail = os.path.split(bam_file)
    chroms = get_contigs_list(arguments['contigs'])
    chroms_without_prefix = [str(i).replace( 'chr', '') for i in chroms]
    for i, seq_elem in enumerate(seq_dict):
        region[i] = seq_elem['LN']
        chrs[i] = seq_elem['SN']
        start=0
        end=bin_size
        if (chrs[i] in chroms) or (chrs[i] in chroms_without_prefix):
            for c in range(0,region[i],bin_size):
                if end > region[i]:
                    bed.append([tail, chrs[i], start, region[i]])
                else:
                    bed.append([tail, chrs[i], start, end])
                start=end+1
                end+=bin_size
    return bed

def chromosomes_sorter(label):
    from itertools import takewhile
    # Strip "chr" prefix
    chrom = (label[3:] if label.lower().startswith('chr')
             else label)
    if chrom in ('X', 'Y'):
        key = (1000, chrom)
    else:
        # Separate numeric and special chromosomes
        nums = ''.join(takewhile(str.isdigit, chrom))
        chars = chrom[len(nums):]
        nums = int(nums) if nums else 0
        if not chars:
            key = (nums, '')
        elif len(chars) == 1:
            key = (2000 + nums, chars)
        else:
            key = (3000 + nums, chars)
    return key

def csv_df_chromosomes_sorter(path, names, sept='\t'):
    dataframe = pd.read_csv(path, sep=sept, names=names)
    dataframe['chr'] = dataframe['chr'].astype(str)
    #if not dataframe['chr'].iloc[0].startswith('chr'):
    #    dataframe['chr'] = 'chr' + dataframe['chr'].astype(str)
    dataframe.sort_values(by=['chr', names[1]], ascending=[True, True], inplace=True)
    return dataframe.reindex(dataframe.chr.apply(chromosomes_sorter).sort_values(kind='mergesort').index)

def df_chromosomes_sorter(dataframe, names, sept='\t'):
    dataframe['chr'] = dataframe['chr'].astype(str)
    #if not dataframe['chr'].iloc[0].startswith('chr'):
    #    dataframe['chr'] = 'chr' + dataframe['chr'].astype(str)
    dataframe.sort_values(by=['chr', names[1]], ascending=[True, True], inplace=True)
    return dataframe.reindex(dataframe.chr.apply(chromosomes_sorter).sort_values(kind='mergesort').index)

def write_segments_coverage(coverage_segments, output):
    with open('data/' + output, 'a') as fp:
        for items in coverage_segments:
            if not items == None:
                fp.write("%s\n" % items)

def seperate_dfs_coverage(arguments, df, haplotype_1_values_updated, haplotype_2_values_updated, unphased):
    if arguments['without_phasing']:
        return df[['chr', 'start', 'end', 'coverage']].copy()
    else:
        df_hp1 = df[['chr', 'start','end', 'hp1']].copy()
        df_hp2 = df[['chr', 'start','end', 'hp2']].copy()
        df_unphased = df[['chr', 'start','end', 'hp3']].copy()
        df_hp1['hp1'] = haplotype_1_values_updated
        df_hp2['hp2'] = haplotype_2_values_updated
        df_unphased['hp3'] = unphased
        return df_hp1, df_hp2, df_unphased


def get_snps_frquncies_coverage_from_bam(df, chrom):
    df = df[df['chr'] == chrom]
    #df = dict(tuple(df.groupby('hp')))
    haplotype_1_position = df.pos.values.tolist()
    haplotype_1_coverage = df.freq_value_b.values.tolist()
    haplotype_2_position = df.pos.values.tolist()
    haplotype_2_coverage = df.freq_value_a.values.tolist()

    return haplotype_1_position, haplotype_1_coverage, haplotype_2_position, haplotype_2_coverage

def detect_alter_loh_regions(arguments, event, chrom, ref_ends, haplotype_1_values, haplotype_2_values, unphased_reads_values, starts, ends):
    if ends and ends[-1] > ref_ends[-1]:
        ends[-1] = ref_ends[-1]

    region_starts = []
    region_ends = []
    #print(starts)
    #print(ends)
    for i, (start,end) in enumerate(zip(starts,ends)):
        if end - start > 2000000:
            region_starts.append(start)
            region_ends.append(end)

    #print(region_starts)
    #print(region_ends)
    if region_starts:
        dict = []
        for i in range(len(region_starts)):
            dict.append((chrom + '\t' + str(region_starts[i]) + '\t' + str(region_ends[i]) + '\t' + event))
        write_segments_coverage(dict, arguments['genome_name'] + '_events_segments.csv')

    for j, (starts,ends) in enumerate(zip(region_starts, region_ends)):
        if mean_values(haplotype_1_values, starts - 1, starts - 4) > mean_values(haplotype_2_values, starts - 1, starts - 4):
            for i in range(starts//arguments['bin_size'],ends//arguments['bin_size']):
                    haplotype_1_values[i] = haplotype_1_values[i] + haplotype_2_values[i] + unphased_reads_values[i]
                    haplotype_2_values[i] = 0
                    unphased_reads_values[i] = 0
        else:
            for i in range(starts // arguments['bin_size'], ends // arguments['bin_size']):
                haplotype_2_values[i] = haplotype_1_values[i] + haplotype_2_values[i] + unphased_reads_values[i]
                haplotype_1_values[i] = 0
                unphased_reads_values[i] = 0

    return haplotype_1_values, haplotype_2_values, unphased_reads_values


def mean_values(selected_list, start_index, end_index):
    result = []
    for i in range(end_index - start_index):
        try:
            result.append(selected_list[start_index + i])
        except IndexError:
            break
    if result:
        return np.mean(result)
    else:
        return 0.0

def infer_missing_phaseblocks(bam_file, chrom, ref_start_values_phasesets, ref_end_values_phasesets):

    head, tail = os.path.split(bam_file)
    bed = []
    for i in range(len(ref_start_values_phasesets)-1):
        if ref_start_values_phasesets[i +1] - ref_end_values_phasesets[i] > 1:
            bed.append([tail, chrom, ref_end_values_phasesets[i] + 1, ref_start_values_phasesets[i+1] - 1])

    return bed

def update_phasesets_coverage_with_missing_phasesets(chroms, csv_df_phasesets, bam, coverage_histograms):
    coverage = []
    for index, chrom in enumerate(chroms):
        csv_df_phaseset = csv_df_phasesets[csv_df_phasesets['chr'] == chrom]
        ref_start_values_phasesets = csv_df_phaseset.start.values.tolist()
        ref_end_values_phasesets = csv_df_phaseset.end.values.tolist()

        missing_ps_segments = infer_missing_phaseblocks(bam, chrom, ref_start_values_phasesets, ref_end_values_phasesets)
        missing_ps_coverage = get_segments_coverage(missing_ps_segments, coverage_histograms)
        coverage.extend(missing_ps_coverage)

    return coverage