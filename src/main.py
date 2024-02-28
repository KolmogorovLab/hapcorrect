#!/usr/bin/env python3

import sys
import shutil

import numpy
import pysam
import argparse
import logging
import os
import pandas as pd

from multiprocessing import Pool
from collections import defaultdict

from process_bam import get_all_reads_parallel, update_coverage_hist, get_segments_coverage, haplotype_update_all_bins_parallel, get_snps_frequencies
from process_vcf import vcf_parse_to_csv_for_het_phased_snps_phasesets, get_snp_frequencies_segments, snps_frequencies_chrom_mean, rephase_vcf
from phase_correction import generate_phasesets_bins, phaseblock_flipping, phase_correction_centers, contiguous_phaseblocks, detect_centromeres, flip_phaseblocks_contigous
from utils import get_chromosomes_bins, write_segments_coverage, csv_df_chromosomes_sorter, get_snps_frquncies_coverage_from_bam, \
                    infer_missing_phaseblocks, df_chromosomes_sorter, is_phasesets_check_simple_heuristics
from extras import get_contigs_list
from plots import plot_coverage_data, change_point_detection, plot_coverage_data_after_correction
from cpd import cpd_positions_means

def main():
    # default tunable parameters
    MAX_READ_ERROR = 0.1
    MIN_MAPQ = 10
    MIN_SV_SIZE = 50
    BIN_SIZE = 50000
    MAX_CUT_THRESHOLD = 100
    MIN_ALIGNED_LENGTH = 5000

    SAMTOOLS_BIN = "samtools"
    BCFTOOLS_BIN = "bcftools"

    DEFAULT_CONTIGS = 'chr1-22' #('chr1-22' '1-22') ('chr1-22,chrX' '1-22,X')

    parser = argparse.ArgumentParser \
        (description="Plot coverage and copy number profiles from a bam and phased VCF files")

    parser.add_argument("--target-bam", dest="target_bam",
                        metavar="path", required=True, default=None, nargs="+",
                        help="path to one or1/4 multiple target haplotagged bam files (must be indexed)")
    parser.add_argument("--control-bam", dest="control_bam",
                        metavar="path", required=False, default=None, nargs="+",
                        help="path to one or multiple control haplotagged bam files (must be indexed)")
    parser.add_argument("--reference", dest="reference",
                        metavar="path", required=False, default=None,
                        help="path to reference")

    parser.add_argument("--out-dir-plots", dest="out_dir_plots",
                        default=None, required=True,
                        metavar="path", help="Output plots directory")

    parser.add_argument("--phased-vcf", dest="phased_vcf",
                        metavar="path", required=True, default=None,
                        help="Path to phased vcf")
    parser.add_argument("--phased-vcf-snps-freqs", dest="phased_vcf_snps_freqs",
                        metavar="path", required=False, default=None,
                        help="Path to phased vcf to plot snps frequencies coverage")

    parser.add_argument("--genome-name", dest="genome_name",
                        required=True, default=None,
                        help="Genome sample/cellline name to be displayed on plots")
    parser.add_argument("--contigs", dest="contigs",
                        required=False, default=DEFAULT_CONTIGS,
                        help="List of contigs (choromosomes) to be included in the plots [e.g., chr1-22,X,Y]")

    parser.add_argument("--bin-size", "--bin_size", dest="bin_size",
                        default=BIN_SIZE, metavar="int", type=int, help="coverage (readdepth) bin size [50k]")
    parser.add_argument("--cut-threshold", "--cut_threshold", dest="cut_threshold",
                        default=MAX_CUT_THRESHOLD, metavar="int", type=int, help="Maximum cut threshold for coverage (readdepth) [100]")

    parser.add_argument("--min-aligned-length", "--min_aligned_length", dest="min_aligned_length",
                        default=MIN_ALIGNED_LENGTH, metavar="int", type=int, help="Minimum aligned reads length [5000]")

    parser.add_argument('--pdf-enable',  dest="pdf_enable", required=False,
                        default=False, help="Enabling PDF output coverage plots")
    parser.add_argument('--html-enable',  dest="html_enable", required=False,
                        default=True, help="Enabling HTML output coverage plots")

    parser.add_argument('--unphased-reads-coverage-enable',  dest="unphased_reads_coverage_enable", required=False,
                        default=False, help="Enabling unphased reads coverage output in plots")

    parser.add_argument('--phaseblock-flipping-enable',  dest="phaseblock_flipping_enable", required=False,
                        default=False, help="Enabling phaseblock flipping in coverage plots")
    parser.add_argument('--phaseblocks-enable',  dest="phaseblocks_enable", required=False,
                        default=True, help="Enabling phaseblocks display in coverage plots")
    parser.add_argument('--het-phased-snps-freq-enable',  dest="het_phased_snps_freq_enable", required=False,
                        default=False, help="Enabling hetrozygous phased snps frequencies in coverage plots")
    parser.add_argument('--snps-freq-vcf-enable',  dest="snps_freq_vcf_enable", required=False,
                        default=False, help="Enabling snps frequencies in coverage plots")

    parser.add_argument("-t", "--threads", dest="threads",
                        default=1, metavar="int", type=int, help="number of parallel threads [8]")
    parser.add_argument('--dryrun', dest="dryrun", required=False,
                        default=False, help="Enabling dryrun")

    parser.add_argument("--max-read-error", dest="max_read_error",
                        default=MAX_READ_ERROR, metavar="float", type=float,
                        help=f"maximum base alignment error [{MAX_READ_ERROR}]")
    parser.add_argument("--min-mapq", dest="min_mapping_quality",
                        default=MIN_MAPQ, metavar="int", type=int,
                        help=f"minimum mapping quality for aligned segment [{MIN_MAPQ}]")

    args = parser.parse_args()

    arguments = {
        "target_bam": args.target_bam,
        "reference": args.reference,
        "phaseblock_flipping_enable": args.phaseblock_flipping_enable,
        "unphased_reads_coverage_enable": args.unphased_reads_coverage_enable,
        "phaseblocks_enable": args.phaseblocks_enable,
        "het_phased_snps_freq_enable": args.het_phased_snps_freq_enable,
        "snps_freq_vcf_enable": args.snps_freq_vcf_enable,
        "out_dir_plots": args.out_dir_plots,
        "phased_vcf": args.phased_vcf,
        "phased_vcf_snps_freqs": args.phased_vcf_snps_freqs,
        "genome_name": args.genome_name,
        "bin_size": args.bin_size,
        "pdf_enable": args.pdf_enable,
        "cut_threshold": args.cut_threshold,
        "contigs": args.contigs,
        "threads": args.threads,
        "dryrun":args.dryrun,
        "min_aligned_length": args.min_aligned_length,
    }
    logging.basicConfig(level=logging.DEBUG)

    if args.control_bam is None:
        args.control_bam = []
    all_bams = args.target_bam + args.control_bam
    target_genomes = set(os.path.basename(b) for b in args.target_bam)
    control_genomes = set(os.path.basename(b) for b in args.control_bam)

    if not shutil.which(SAMTOOLS_BIN):
        print("samtools not found", file=sys.stderr)
        return 1

    # TODO: check that all bams have the same reference
    first_bam = all_bams[0]
    ref_lengths = None
    with pysam.AlignmentFile(first_bam, "rb") as a:
        ref_lengths = dict(zip(a.references, a.lengths))

    if not os.path.isdir(args.out_dir_plots):
        os.mkdir(args.out_dir_plots)
    #if not os.path.isdir('data'):
    if os.path.exists('data'):
        shutil.rmtree('data')
        os.mkdir('data')
    else:
        os.mkdir('data')

    thread_pool = Pool(args.threads)

    segments_by_read = defaultdict(list)
    genome_ids = []
    for bam_file in all_bams:
        genome_id = os.path.basename(bam_file)
        genome_ids.append(genome_id)
        print("Parsing reads from", genome_id, file=sys.stderr)
        segments_by_read_bam = get_all_reads_parallel(bam_file, thread_pool, ref_lengths,
                                                      args.min_mapping_quality, genome_id, MIN_SV_SIZE)
        segments_by_read.update(segments_by_read_bam)
        print("Parsed {0} segments".format(len(segments_by_read_bam)), file=sys.stderr)

    logging.info('Computing coverage histogram')
    coverage_histograms = update_coverage_hist(genome_ids, ref_lengths, segments_by_read, args.min_mapping_quality,
                                               args.max_read_error, arguments)
    del segments_by_read
    chroms = get_contigs_list(arguments['contigs'])

    if arguments['dryrun']:
        csv_df_phasesets = csv_df_chromosomes_sorter('/home/rezkuh/gits/data/' + arguments['genome_name'] + '/coverage_ps.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
        #csv_df_phasesets_missing = csv_df_chromosomes_sorter('/home/rezkuh/gits/data/' + arguments['genome_name'] + '/coverage_ps_missing.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
        #csv_df_phasesets = df_chromosomes_sorter(pd.concat([csv_df_phasesets, csv_df_phasesets_missing], ignore_index=True), ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])

        csv_df_coverage = csv_df_chromosomes_sorter('/home/rezkuh/gits/data/' + arguments['genome_name'] + '/coverage.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
    else:
        logging.info('Computing coverage for bins')
        segments = get_chromosomes_bins(args.target_bam[0], arguments['bin_size'], arguments)
        segments_coverage = get_segments_coverage(segments, coverage_histograms)
        logging.info('Writing coverage for bins')
        write_segments_coverage(segments_coverage, 'coverage.csv')

        logging.info('Parsing phaseblocks information')
        output_phasesets_file_path = vcf_parse_to_csv_for_het_phased_snps_phasesets(arguments['phased_vcf'])
        phasesets_segments = generate_phasesets_bins(args.target_bam[0], output_phasesets_file_path, arguments['bin_size'], arguments) #TODO update for multiple bam files
        logging.info('Computing coverage for phaseblocks')
        phasesets_coverage = get_segments_coverage(phasesets_segments, coverage_histograms)

        logging.info('Writing coverage for phaseblocks')
        write_segments_coverage(phasesets_coverage, 'coverage_ps.csv')

        logging.info('Loading coverage (bins) and coverage (phaseblocks) files...')
        csv_df_phasesets = csv_df_chromosomes_sorter('data/coverage_ps.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])
        csv_df_coverage = csv_df_chromosomes_sorter('data/coverage.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])

        #Missing phaseblocks and coverage added
        #phasesets_coverage_missing = update_phasesets_coverage_with_missing_phasesets(chroms, csv_df_phasesets, args.target_bam[0], coverage_histograms)
        #write_segments_coverage(phasesets_coverage_missing, 'coverage_ps_missing.csv')
        #csv_df_phasesets_missing = csv_df_chromosomes_sorter('data/coverage_ps_missing.csv', ['chr', 'start', 'end', 'hp1', 'hp2', 'hp3'])

        del coverage_histograms

    get_snp_frequencies_segments(arguments, arguments['target_bam'][0], thread_pool)
    df_snps_frequencies = csv_df_chromosomes_sorter('data/snps_frequencies.csv', ['chr', 'pos', 'freq_value_a', 'hp_a', 'freq_value_b', 'hp_b'])
    df_snps_frequencies = df_snps_frequencies.drop(df_snps_frequencies[(df_snps_frequencies.chr == "chrX") | (df_snps_frequencies.chr == "chrY")].index)

    filename = f"{os.path.join(arguments['out_dir_plots'], 'COVERAGE.html')}"
    html_graphs = open(filename, 'w')
    html_graphs.write("<html><head></head><body>" + "\n")

    for index, chrom in enumerate(chroms):

        if chrom in chroms:# and (chrom == '1' or chrom == '18'):# and (chrom == 'chr5' or chrom == 'chr16'):
            logging.info('Loading coverage (bins) and coverage (phaseblocks) datasets for ' + chrom)
            csv_df_phaseset = csv_df_phasesets[csv_df_phasesets['chr'] == chrom]
            haplotype_1_values_phasesets = csv_df_phaseset.hp1.values.tolist()
            haplotype_2_values_phasesets = csv_df_phaseset.hp2.values.tolist()
            ref_start_values_phasesets = csv_df_phaseset.start.values.tolist()
            ref_end_values_phasesets = csv_df_phaseset.end.values.tolist()

            csv_df_coverage_chrom = csv_df_coverage[csv_df_coverage['chr'] == chrom]
            unphased_reads_values = csv_df_coverage_chrom.hp3.values.tolist()
            haplotype_1_values = csv_df_coverage_chrom.hp1.values.tolist()
            haplotype_2_values = csv_df_coverage_chrom.hp2.values.tolist()
            ref_start_values = csv_df_coverage_chrom.start.values.tolist()
            ref_end_values = csv_df_coverage_chrom.end.values.tolist()

            snps_haplotype1_mean, snps_haplotype2_mean  = snps_frequencies_chrom_mean(df_snps_frequencies, ref_start_values, chrom)
            plot_coverage_data(html_graphs, arguments, chrom, ref_start_values, ref_end_values, snps_haplotype1_mean,
                               snps_haplotype2_mean, unphased_reads_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets, "without_phase_correction")
            ##################################
            # change_point_detection(snps_haplotype1_mean, ref_start_values, ref_end_values, arguments, chrom,
            #                        html_graphs, 1, color='#6A5ACD')
            # change_point_detection(snps_haplotype2_mean, ref_start_values, ref_end_values, arguments, chrom,
            #                        html_graphs, 2, color='#2E8B57')

            # cpd_haplotype1_mean, cpd_haplotype1_start, cpd_haplotype1_end,  cpd_haplotype2_mean, cpd_haplotype2_start, cpd_haplotype2_end = cpd_positions_means(snps_haplotype1_mean, snps_haplotype2_mean, arguments)

            # change_point_detection(snps_haplotype1_mean, ref_start_values, ref_end_values, arguments, chrom,
            #                        html_graphs, 1, color='#6A5ACD')
            # change_point_detection(snps_haplotype2_mean, ref_start_values, ref_end_values, arguments, chrom,
            #                        html_graphs, 2, color='#2E8B57')

            # phase_correction_centers(arguments, cpd_haplotype1_mean, cpd_haplotype1_start, cpd_haplotype1_end,  cpd_haplotype2_mean, cpd_haplotype2_start, cpd_haplotype2_end, snps_haplotype1_mean, snps_haplotype2_mean)
            # print(cpd_haplotype1_mean)
            # print([(cpd_haplotype1_end - cpd_haplotype1_start) //50000 for cpd_haplotype1_start, cpd_haplotype1_end in zip(cpd_haplotype1_start, cpd_haplotype1_end)])
            # print(cpd_haplotype2_mean)
            # print([(cpd_haplotype2_end - cpd_haplotype2_start) // 50000 for cpd_haplotype2_start, cpd_haplotype2_end in
            #       zip(cpd_haplotype2_start, cpd_haplotype2_end)])

            # print(cpd_haplotype2_mean,cpd_haplotype2_end-cpd_haplotype2_start)

            ##################################
            #phaseblocks flipping main algo
            #is_phasesets_check_simple_heuristics(ref_start_values_phasesets, ref_end_values_phasesets)
            snps_haplotype1_mean, snps_haplotype2_mean, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets = \
                phaseblock_flipping(chrom, arguments, False, snps_haplotype1_mean, snps_haplotype2_mean, ref_start_values, ref_start_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)

            # #detect centromeres
            ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets = detect_centromeres(ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets, snps_haplotype1_mean, snps_haplotype2_mean, ref_start_values, arguments['bin_size'])
            #infer missing phaseblocks
            ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets = infer_missing_phaseblocks(ref_start_values, ref_end_values, ref_start_values_phasesets, ref_end_values_phasesets, haplotype_1_values_phasesets, haplotype_2_values_phasesets, snps_haplotype1_mean, snps_haplotype2_mean, arguments['bin_size'])
            #create more contiguous phaseblocks
            haplotype_1_values_phasesets_conti, ref_start_values_phasesets_hp1, ref_end_values_phasesets_hp1, haplotype_2_values_phasesets_conti, ref_start_values_phasesets_hp2, ref_end_values_phasesets_hp2 = contiguous_phaseblocks(
                haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)

            #flip phaseblocks based on more contiguous phaseblocks
            haplotype_1_values_phasesets, haplotype_2_values_phasesets, snps_haplotype1_mean, snps_haplotype2_mean = flip_phaseblocks_contigous(chrom, arguments,
                haplotype_1_values_phasesets_conti, ref_start_values_phasesets_hp1, ref_end_values_phasesets_hp1,
                haplotype_2_values_phasesets_conti, ref_start_values_phasesets_hp2, ref_end_values_phasesets_hp2,
                ref_start_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets,
                ref_start_values_phasesets, ref_end_values_phasesets, snps_haplotype1_mean, snps_haplotype2_mean)

            # #plot resultant
            plot_coverage_data(html_graphs, arguments, chrom, ref_start_values, ref_end_values, snps_haplotype1_mean,
                               snps_haplotype2_mean, unphased_reads_values, haplotype_1_values_phasesets,
                               haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets,
                               "with_phase_correction_0")
            ##################################

            ##################################
            haplotype_1_values_phasesets_conti, ref_start_values_phasesets_hp1, ref_end_values_phasesets_hp1, haplotype_2_values_phasesets_conti, ref_start_values_phasesets_hp2, ref_end_values_phasesets_hp2 = contiguous_phaseblocks(
                haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets)
            plot_coverage_data_after_correction(html_graphs, arguments, chrom, ref_start_values, ref_end_values, snps_haplotype1_mean,
                               snps_haplotype2_mean, unphased_reads_values, haplotype_1_values_phasesets_conti, haplotype_2_values_phasesets_conti, ref_start_values_phasesets_hp1, ref_end_values_phasesets_hp1, ref_start_values_phasesets_hp2, ref_end_values_phasesets_hp2, "phase_correction")
    html_graphs.write("</body></html>")

    # TODO call edit VCF functionality from process_vcf() here
    csv_df_phase_change_segments = csv_df_chromosomes_sorter('data/' + arguments['genome_name'] + '_phase_change_segments.csv', ['chr', 'start', 'end'])
    logging.info('VCF edit for phase change segments')
    out_vcf = os.path.join(arguments['out_dir_plots'], 'rephased_vcf.vcf.gz')
    rephase_vcf(csv_df_phase_change_segments, arguments["phased_vcf"],out_vcf)


if __name__ == "__main__":
    main()
