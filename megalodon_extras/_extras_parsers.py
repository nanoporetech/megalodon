import sys
import argparse

from megalodon import megalodon_helper as mh


####################
# Aggregate Parser #
####################

def get_parser_aggregate_run():
    parser = argparse.ArgumentParser(
        'Aggregate per-read, per-site statistics from previous megalodon ' +
        'call.')

    out_grp = parser.add_argument_group('Output Arguments')
    out_grp.add_argument(
        '--outputs', nargs='+',
        default=[mh.VAR_NAME, mh.MOD_NAME],
        choices=[mh.VAR_NAME, mh.MOD_NAME],
        help='Output type(s) to produce. Default: %(default)s')
    out_grp.add_argument(
        '--megalodon-directory', default='megalodon_results',
        help='Megalodon output directory containing per-read database(s) ' +
        'where aggregated results will be added. Default: %(default)s')
    out_grp.add_argument(
        '--output-suffix', default='re_aggregated',
        help='Suffix to apply to aggregated results, to avoid ' +
        'overwriting results. Default: %(default)s')
    out_grp.add_argument(
        '--read-ids-filename',
        help='File containing read ids to process (one per ' +
        'line). Default: All reads')

    var_grp = parser.add_argument_group('Sequence Variant Arguments')
    var_grp.add_argument(
        '--haploid', action='store_true',
        help='Compute sequence variant aggregation for haploid genotypes. ' +
        'Default: diploid')
    var_grp.add_argument(
        '--heterozygous-factors', type=float, nargs=2,
        default=[mh.DEFAULT_SNV_HET_FACTOR, mh.DEFAULT_INDEL_HET_FACTOR],
        help='Bayesian prior factor for snv and indel heterozygous calls ' +
        '(compared to 1.0 for hom ref/alt). Default: %(default)s')
    var_grp.add_argument(
        '--write-vcf-log-probs', action='store_true',
        help='Write alt log prbabilities out in non-standard VCF field.')

    mod_grp = parser.add_argument_group('Modified Base Arguments')
    mod_grp.add_argument(
        '--mod-aggregate-method', choices=list(mh.MOD_AGG_METHOD_NAMES),
        default=mh.MOD_BIN_THRESH_NAME,
        help='Modified base aggregation method. Default: %(default)s')
    mod_grp.add_argument(
        '--mod-binary-threshold', type=float, nargs=1,
        default=mh.DEFAULT_MOD_BINARY_THRESH,
        help='Threshold for modified base aggregation (probability of ' +
        'modified/canonical base). Only applicable for ' +
        '"--mod-aggregate-method binary_threshold". Default: %(default)s')
    mod_grp.add_argument(
        '--mod-output-formats', nargs='+',
        default=[mh.MOD_BEDMETHYL_NAME, ],
        choices=tuple(mh.MOD_OUTPUT_FMTS.keys()),
        help='Modified base aggregated output format(s). Default: %(default)s')
    mod_grp.add_argument(
        '--write-mod-log-probs', action='store_true',
        help='Write per-read modified base log probabilities ' +
        'out in non-standard modVCF field.')

    misc_grp = parser.add_argument_group('Miscellaneous Arguments')
    misc_grp.add_argument(
        '--processes', type=int, default=1,
        help='Number of parallel processes. Default: %(default)d')
    misc_grp.add_argument(
        '--suppress-progress', action='store_true',
        help='Suppress progress bar output.')

    return parser


#####################
# Calibrate Parsers #
#####################

def get_parser_calibrate_modified_bases():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--ground-truth-llrs', default='mod_calibration_statistics.npz',
        help='Ground truth log-likelihood ratio statistics (produced by ' +
        '`megalodon_extras calibrate generate_modified_base_stats`). ' +
        'Default: %(default)s')
    parser.add_argument(
        '--max-input-llr', type=int, default=mh.DEFAULT_CALIB_SMOOTH_MAX,
        help='Maximum log-likelihood ratio to compute calibration. ' +
        'Default: %(default)d')
    parser.add_argument(
        '--num-calibration-values', type=int,
        default=mh.DEFAULT_CALIB_SMOOTH_NVALS,
        help='Number of discrete calibration values to compute. ' +
        'Default: %(default)d')
    parser.add_argument(
        '--smooth-bandwidth', type=float,
        default=mh.DEFAULT_CALIB_SMOOTH_BW,
        help='Smoothing bandwidth. Default: %(default)f')
    parser.add_argument(
        '--min-density', type=float, default=mh.DEFAULT_CALIB_MIN_DENSITY,
        help='Minimum density value to compute calibration. This value ' +
        'dynamically adjusts [--max-input-llr] when it is too large. ' +
        'Default: %(default)f')
    parser.add_argument(
        '--diff-epsilon', type=float, default=mh.DEFAULT_CALIB_DIFF_EPS,
        help='Epsilon to determine when the likelihood ratio has plateaued. ' +
        'Default: %(default)f')
    parser.add_argument(
        '--llr-clip-buffer', type=int,
        default=mh.DEFAULT_CALIB_LLR_CLIP_BUFFER,
        help='Clipped buffer when determining range for computed ' +
        'calibration log likelihood ratios. Default: %(default)d')
    parser.add_argument(
        '--out-filename', default='megalodon_mod_calibration.npz',
        help='Filename to output calibration values. Default: %(default)s')
    parser.add_argument(
        '--out-pdf',
        help='Output pdf filename for modified base calibration ' +
        'visualization. Default: Do not produce plot.')
    parser.add_argument(
        '--pdf-prob-thresholds', nargs=3, type=float,
        default=[0.75, 0.8, 0.85],
        help='Probability thresholds to mark on output pdf. ' +
        'Default: %(default)s')
    parser.add_argument(
        '--plot-without-prob-thresholds', action='store_true',
        help='Do not include probability thresholds in plot(s).')
    parser.add_argument(
        '--processes', type=int, default=1,
        help='Number of processing cores to use for density smoothing ' +
        'computation. Default: %(default)d')
    parser.add_argument(
        '--overwrite', action='store_true',
        help='Overwrite --out-filename if it exists.')

    return parser


def get_parser_calibrate_merge_modified_bases():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'modified_base_calibration_files', nargs='+', metavar='MOD_CALIB_FN',
        help='Modified base calibration filenames. For modified bases ' +
        'included in more than one file the values from the first file ' +
        'listed will be used.')
    parser.add_argument(
        '--out-filename', default='megalodon_mod_calibration.npz',
        help='Filename to output calibration values. Default: %(default)s')
    parser.add_argument(
        '--overwrite', action='store_true',
        help='Overwrite --out-filename if it exists.')

    return parser


def get_parser_calibrate_variants():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--ground-truth-llrs', default='variant_calibration_statistics.txt',
        help='Ground truth log-likelihood ratio statistics (produced by ' +
        '`megalodon_extras calibrate generate_variant_stats`). ' +
        'Default: %(default)s')
    parser.add_argument(
        '--max-input-llr', type=int, default=mh.DEFAULT_CALIB_SMOOTH_MAX,
        help='Maximum log-likelihood ratio to compute calibration. ' +
        'Default: %(default)d')
    parser.add_argument(
        '--num-calibration-values', type=int,
        default=mh.DEFAULT_CALIB_SMOOTH_NVALS,
        help='Number of discrete calibration values to compute. ' +
        'Default: %(default)d')
    parser.add_argument(
        '--smooth-bandwidth', type=float,
        default=mh.DEFAULT_CALIB_SMOOTH_BW,
        help='Smoothing bandwidth. Default: %(default)f')
    parser.add_argument(
        '--min-density', type=float, default=mh.DEFAULT_CALIB_MIN_DENSITY,
        help='Minimum density value to compute calibration. This value ' +
        'dynamically adjusts [--max-input-llr] when it is too large. ' +
        'Default: %(default)f')
    parser.add_argument(
        '--diff-epsilon', type=float, default=mh.DEFAULT_CALIB_DIFF_EPS,
        help='Epsilon to determine when the likelihood ratio has plateaued. ' +
        'Default: %(default)f')
    parser.add_argument(
        '--llr-clip-buffer', type=int,
        default=mh.DEFAULT_CALIB_LLR_CLIP_BUFFER,
        help='Clipped buffer when determining range for computed ' +
        'calibration log likelihood ratios. Default: %(default)d')
    parser.add_argument(
        '--processes', type=int, default=1,
        help='Number of processing cores to use for density smoothing ' +
        'computation. Default: %(default)d')
    parser.add_argument(
        '--out-filename', default='megalodon_variant_calibration.npz',
        help='Filename to output calibration values. Default: %(default)s')
    parser.add_argument(
        '--out-pdf',
        help='Output pdf filename for modified base calibration ' +
        'visualization. Default: Do not produce plot.')
    parser.add_argument(
        '--overwrite', action='store_true',
        help='Overwrite --out-filename if it exists.')

    return parser


def get_parser_calibrate_generate_modified_bases_stats():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'megalodon_results_dir',
        help='Output directory from Megalodon with mappings and ' +
        'per_read_mods in outputs.')
    parser.add_argument(
        '--control-megalodon-results-dir',
        help='Megalodon output directory with modified base control sample.')
    parser.add_argument(
        '--exclude-modified-bases', nargs='+',
        help='Set of modified bases (single letter codes) to exclude.')
    parser.add_argument(
        '--modified-bases-set', nargs='+',
        help='Only process these modified bases (single letter codes).')
    parser.add_argument(
        '--ground-truth-data',
        help='Ground truth csv with (chrm, strand, pos, is_mod) values.')
    parser.add_argument(
        '--strand-specific-sites', action='store_true',
        help='Sites in --ground-truth-data are strand-specific. If not ' +
        'set, strand is ignored.')
    parser.add_argument(
        '--out-filename', default='mod_calibration_statistics.npz',
        help='Output filename for text summary. Default: %(default)s')
    parser.add_argument(
        '--quiet', action='store_true',
        help='Suppress progress information.')

    return parser


def get_parser_calibrate_generate_variants_stats():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'fast5s_dir',
        help='Directory containing raw fast5 (will be searched recursively).')

    pyg_grp = parser.add_argument_group('Guppy Backend Arguments')
    pyg_grp.add_argument(
        '--guppy-config', default=mh.DEFAULT_GUPPY_CFG,
        help='Guppy config. Default: %(default)s')
    pyg_grp.add_argument(
        '--guppy-server-path', default=mh.DEFAULT_GUPPY_SERVER_PATH,
        help='Path to guppy server executable. Default: %(default)s')
    pyg_grp.add_argument(
        '--guppy-server-port', type=int,
        help='Guppy server port. Default: Guppy auto')
    pyg_grp.add_argument(
        '--guppy-params',
        help='Extra guppy server parameters. Main purpose for optimal ' +
        'performance based on compute environment. Quote parameters to ' +
        'avoid them being parsed by megalodon.')
    pyg_grp.add_argument(
        '--guppy-timeout', type=float, default=mh.DEFAULT_GUPPY_TIMEOUT,
        help='Timeout to wait for guppy server to call a single read in ' +
        'seconds. Default: %(default)f')
    pyg_grp.add_argument(
        '--guppy-logs-output-directory', default='guppy_logs',
        help='Directory to output guppy logs. Default: %(default)s')

    map_grp = parser.add_argument_group('Mapping Arguments')
    map_grp.add_argument(
        '--reference',
        help='Reference FASTA file used for mapping called reads.')

    out_grp = parser.add_argument_group('Output Arguments')
    out_grp.add_argument(
        '--output', default='variant_calibration_statistics.txt',
        help='Filename to output statistics. Default: %(default)s')
    out_grp.add_argument(
        '--num-reads', type=int,
        help='Number of reads to process. Default: All reads')
    out_grp.add_argument(
        '--read-ids-filename',
        help='File containing read ids to process (one per ' +
        'line). Default: All reads')

    misc_grp = parser.add_argument_group('Miscellaneous Arguments')
    misc_grp.add_argument(
        '--devices', nargs='+',
        help='GPU devices for guppy or taiyaki basecalling backends.')
    misc_grp.add_argument(
        '--not-recursive', action='store_true',
        help='Only search for fast5 read files directly found within the ' +
        'fast5 directory. Default: search recursively')
    misc_grp.add_argument(
        '--processes', type=int, default=1,
        help='Number of parallel processes. Default: %(default)d')
    misc_grp.add_argument(
        '--suppress-progress', action='store_true',
        help='Suppress progress bar.')
    misc_grp.add_argument(
        '--compute-false-reference-scores', action='store_true',
        help='Compute scores given a false reference. Default: compute ' +
        'all scores with ground truth correct reference.' +
        '***** Experimental feature, may contain bugs *****.')

    return parser


def get_parser_merge_modified_bases():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'megalodon_results_dirs', nargs='+',
        help='Megalodon directories with per_read_mods in output.')
    parser.add_argument(
        '--output-megalodon-results-dir',
        default='megalodon_merge_mods_results',
        help='Output directory. Cannot exist before this command. ' +
        'Default: %(default)s')
    parser.add_argument(
        '--data-batch-size', type=int, default=100000,
        help='Batch size to insert position and statistics data. ' +
        'Default: %(default)d')
    parser.add_argument(
        '--max-processes', type=int, default=4,
        help='Maximum number of processes to open for reading statistics. ' +
        'Each process must load all output database in memory indices, and ' +
        'thus may incur high memory usage. Default: %(default)d')
    parser.add_argument(
        '--single-process', action='store_true',
        help='Do not use multiprocessing with one input database per process.')
    parser.add_argument(
        '--database-safety', type=int, default=0,
        help='Setting for database performance versus corruption ' +
        'protection. Options: 0 (DB corruption on application crash), ' +
        '1 (DB corruption on system crash), 2 (DB safe mode). ' +
        'Default: %(default)d')
    parser.add_argument(
        '--overwrite', action='store_true',
        help='Overwrite output directory if it exists.')

    return parser


def get_parser_merge_aggregated_modified_bases():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'bed_methyl_files', nargs='+',
        help='Input bedmethyl format files.')
    parser.add_argument(
        '--output-bed-methyl-file',
        default='merged_modified_bases.bed',
        help='Output bedmethyl filename. Cannot exist before this command. ' +
        'Default: %(default)s')
    parser.add_argument(
        '--sorted-inputs', action='store_true',
        help='If input bedmethyl files are sorted, files will be merged ' +
        'without reading full file into memory. Sort order should be ' +
        '`sort -k1,1V -k2,2n`.')

    return parser


def get_parser_merge_variants():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'megalodon_results_dirs', nargs='+',
        help='Output megalodon directories with per_read_vars in output.')
    parser.add_argument(
        '--output-megalodon-results-dir',
        default='megalodon_merge_vars_results',
        help='Output directory. Default: %(default)s')
    parser.add_argument(
        '--overwrite', action='store_true',
        help='Overwrite output directory if it exists.')
    parser.add_argument(
        '--var-locations-on-disk', action='store_true',
        help='Force sequnece variant locations to be stored only within on ' +
        'disk database table. This option will reduce the RAM memory ' +
        'requirement, but may slow processing. Default: ' +
        'Store positions in memory.')

    return parser


#########################
# Modified Base Parsers #
#########################

def get_parser_modified_bases_describe_alphabet():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--log-directory', default='.',
        help='Directory to output megalodon log. Default: current ' +
        'working directory.')

    pyg_grp = parser.add_argument_group('Guppy Backend Arguments')
    pyg_grp.add_argument(
        '--guppy-config', default=mh.DEFAULT_GUPPY_CFG,
        help='Guppy config. Default: %(default)s')
    pyg_grp.add_argument(
        '--guppy-server-path', default=mh.DEFAULT_GUPPY_SERVER_PATH,
        help='Path to guppy server executable. Default: %(default)s')
    pyg_grp.add_argument(
        '--guppy-logs-output-directory', default='guppy_logs',
        help='Directory to output guppy logs. Default: %(default)s')
    pyg_grp.add_argument(
        '--do-not-use-guppy-server', action='store_true',
        help='Use alternative basecalling backend. Either FAST5 ' +
        '(default; requires --post_out when running guppy) or taiyaki ' +
        '(set `--taiyaki-model-filename` to use taiyaki backend).')
    pyg_grp.add_argument(
        '--guppy-params',
        help='Extra guppy server parameters. Main purpose for ' +
        'optimal performance based on compute environment. ' +
        'Quote parameters to avoid them being parsed by megalodon.')

    f5_grp = parser.add_argument_group('FAST5 Backend Arguments')
    f5_grp.add_argument(
        '--fast5s_dir',
        help='Directory containing raw fast5.')

    tai_grp = parser.add_argument_group('Taiyaki Backend Arguments')
    tai_grp.add_argument(
        '--taiyaki-model-filename',
        help='Taiyaki basecalling model checkpoint file.')

    return parser


def get_parser_modified_bases_estimate_threshold():
    parser = argparse.ArgumentParser(
        description='Estimate ideal threshold for marking up modified bases.')

    parser.add_argument(
        'megalodon_results_dir',
        help='Output directory from megalodon with per_read_mods in output.')
    parser.add_argument(
        'mod_base', help='Single letter code for the modified base.')

    parser.add_argument(
        '--fraction-modified', type=float,
        help='Specify fraction of modified calls. Default: Use ' +
        '--mod-percentile most extreme scores to estimate the fraction.')
    parser.add_argument(
        '--mod-percentile', type=float, default=8.0,
        help='Percentile of extreme scores to determine fraction of ' +
        'modified bases. Default: %(default)d')
    parser.add_argument(
        '--num-statistics', type=int,
        help='Number of per-read statistics to use in estimation. ' +
        'Default: All statistics')

    return parser


def get_parser_modified_bases_update_database():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'old_db',
        help='Megalodon version 0.1 modified base data base.')
    parser.add_argument(
        '--new-db', default='megalodon_mods.db',
        help='Output data base name. Should replace ' +
        'per_read_modified_base_calls.db in megalodon results directory in ' +
        'order to process further. Default: %(default)s')

    return parser


def get_parser_modified_bases_split_calls_by_motif():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'reference',
        help='Reference FASTA file. Must include index file ending in fai.')
    parser.add_argument(
        '--motif', nargs=2, action='append', required=True,
        metavar=['MOTIF', 'REL_POS'],
        help='Motif description. Motifs include two values specifying the ' +
        'sequence motif (may include ambiguity codes) and the relative ' +
        'modified position. Multiple `--motif` values should be provided.')
    parser.add_argument(
        '--megalodon-directory', default='megalodon_results',
        help='Megalodon output directory containing per-read modified base ' +
        'database to be split. Default: %(default)s')
    parser.add_argument(
        '--output-suffix', default='split_by_motif',
        help='Suffix to apply to log (stored in input directory). ' +
        'Default: %(default)s')
    parser.add_argument(
        '--output-prefix', default='megalodon_results.split_by_motif',
        help='Prefix for output directories. One directory will be created ' +
        'for each motif with names [--output-prefix].[--motif]. ' +
        'Default: %(default)s')

    return parser


def get_parser_modified_bases_create_ground_truth():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--bed-methyl-files', nargs='+', required=True,
        help='Bed methyl files.')
    parser.add_argument(
        '--coverage-threshold', type=int, default=1,
        help='Only include sites with sufficient coverage. ' +
        'Default: 1 (= All sites)')
    parser.add_argument(
        '--pct-mod-thresholds', type=float, nargs=2, default=[10.0, 90.0],
        help='Lower and upper percent modified thresholds for ground truth ' +
        'modified positions. Default: %(default)s')
    parser.add_argument(
        '--out-csv', default='ground_truth_modifications.csv',
        help='Output filename for ground truth calls. Default: %(default)s')
    parser.add_argument(
        '--strand-offset', type=int,
        help='Offset to combine stranded results. Positive value indicates ' +
        'reverse strand sites have higher position values. Default treat ' +
        'strands independently.')

    return parser


def get_parser_modified_bases_create_motif_bed():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'reference',
        help='Reference FASTA file. Must include index file ending in fai.')
    parser.add_argument(
        '--motif', nargs=2, action='append', required=True,
        metavar=['MOTIF', 'REL_POS'],
        help='Motif description. Motifs include two values specifying the ' +
        'sequence motif (may include ambiguity codes) and the relative ' +
        'modified position. Multiple `--motif` values should be provided.')
    parser.add_argument(
        '--out-filename', default='motif_sites.bed',
        help='Output BED filename. Default: %(default)s')

    return parser


def get_parser_modified_bases_per_site_thresholds():
    parser = argparse.ArgumentParser(
        description='Extract Megalodon modified base score thresholds at ' +
        'each covered site for marking up signal mapping sequences.')

    parser.add_argument(
        'megalodon_results_dir',
        help='Output directory from megalodon with per_read_mods in output.')
    parser.add_argument(
        'ground_truth_bed',
        help='BEDmethyl file containing ground truth fraction modified. ' +
        'File must be sorted (`sort -k1V -k2n`).')

    parser.add_argument(
        '--ground-truth-cov-min', type=int, default=15,
        help='Minimum coverage (both strands) to include a site from ' +
        'ground truth data. Default: %(default)d')
    parser.add_argument(
        '--nanopore-cov-min', type=int, default=30,
        help='Minimum coverage (single strand) to include a site from ' +
        'nanopore data. Default: %(default)d')
    parser.add_argument(
        '--mod-bases', default='m',
        help='Single letter codes for the modified base. For ' +
        'mulitple alternative bases supply all single letter codes ' +
        'with no spaces. Default: %(default)s')
    parser.add_argument(
        '--strand-offset', type=int,
        help='Offset to combine stranded results. Positive value indicates ' +
        'reverse strand sites have higher position values. Default treat ' +
        'strands independently.')
    parser.add_argument(
        '--valid-sites', nargs='+',
        help='BED files containing sites over which to restrict ' +
        'modified base results. Useful when processing full results using a ' +
        'subset of the ground truth (e.g. CG and CH processing). Must be ' +
        'sorted in same order as [ground_truth_bed] (`sort -k1V -k2n`)')
    parser.add_argument(
        '--out-low-coverage-sites', default='low_coverage_sites.bed',
        help='Output filename for sites with low ground truth or nanopore ' +
        'coverage. Default: %(default)s')
    parser.add_argument(
        '--out-per-site-mod-thresholds', default='site_mod_thresholds.bed',
        help='Output filename for per-site megalodon mod scoring ' +
        'thresholds. Default: %(default)s')
    parser.add_argument(
        '--log-filename', default='per_site_thresholds.log',
        help='Output filename for logging. Default: %(default)s')

    parser.add_argument(
        '--batch-size', type=int, default=mh.DEFAULT_BEDMETHYL_BATCH,
        help='Number of sites to include in each batch for processing. ' +
        'Default: %(default)d')
    parser.add_argument(
        '--processes', type=int, default=1,
        help='Number of processes. Default: %(default)d')

    return parser


def get_parser_modified_bases_index_database():
    parser = argparse.ArgumentParser(
        description='Create per-read modified bases calls database index. ' +
        'Can rescue results from unexpected program crashes.')

    parser.add_argument(
        '--megalodon-directory', default='megalodon_results',
        help='Megalodon output directory containing per-read modified bases ' +
        'database to be indexed. Default: %(default)s')
    parser.add_argument(
        '--output-suffix', default='mod_index_database',
        help='Log file suffix. Default: %(default)s')

    return parser


#########################
# Phase Variant Parsers #
#########################

def get_parser_phase_variants_whatshap_filter():
    parser = argparse.ArgumentParser(
        description='Remove variants incompatible with whatshap')
    parser.add_argument(
        'in_vcf', help='Megalodon VCF file')
    parser.add_argument(
        'out_vcf', help='Output VCF file')
    parser.add_argument(
        '--filtered-records', help='File to output filtered records.')

    return parser


def get_parser_phase_variants_extract_haplotype_reads():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'alignment_filename', help='Alignment filename.')
    parser.add_argument(
        'out_basename', help='Basename for read ids output.')

    return parser


def get_parser_phase_variants_merge_haploid_variants():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'diploid_called_variants',
        help='Phased variants from which the diploid calls are derived.')
    parser.add_argument(
        'haplotype1_variants',
        help='Variant file for haplotype 1.')
    parser.add_argument(
        'haplotype2_variants',
        help='Variant file for haplotype 2.')
    parser.add_argument(
        '--out-vcf', default='merged_haploid_variants.vcf',
        help='Output name for VCF. Default: %(default)s')
    parser.add_argument(
        '--force-invalid-variant-processing', action='store_true',
        help='Force processing of mismatching varints. This script is ' +
        'intended only to process variant files produced from the same set ' +
        'of megalodon per-read variant calls. Behavior when processing ' +
        'mismatched variants is not defined.')

    return parser


#######################
# Text Output Parsers #
#######################

def get_parser_per_read_text_modified_bases():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'megalodon_results_dir',
        help='Output directory from megalodon with per_read_mods in output.')
    parser.add_argument(
        '--out-filename',
        help='Output filename for text summary. Default: output into ' +
        'megalodon results directory')

    return parser


def get_parser_per_read_text_variants():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'megalodon_results_dir',
        help='Output directory from megalodon with per_read_variants ' +
        'in output.')
    parser.add_argument(
        '--out-filename',
        help='Output filename for text summary. Default: output into ' +
        'megalodon results directory')

    return parser


####################
# Validate Parsers #
####################

def get_parser_validate_results():
    parser = argparse.ArgumentParser(
        description='Produce per-read results report for sequence mappings ' +
        'and modified bases. Modified base validation requires a ground ' +
        'truth in the form of either a control sample or ground truth ' +
        'modified and unmodified sites within a sample.')

    parser.add_argument(
        'megalodon_results_dirs', nargs='+',
        help='Output directories from megalodon with mappings and ' +
        'optionally per_read_mods in outputs.')

    mod_grp = parser.add_argument_group('Modified Base Arguments')
    mod_grp.add_argument(
        '--control-megalodon-results-dirs', nargs='+',
        help='Megalodon output directories for modified base control ' +
        'sample(s). Could be a PCR or IVT sample. Either a single control ' +
        'for all modified samples or one control sample for each modified ' +
        'sample should be provided.')
    mod_grp.add_argument(
        '--ground-truth-data',
        help='Ground truth csv with (chrm, strand, pos, is_mod) values.')
    mod_grp.add_argument(
        '--valid-sites', nargs=2, action='append',
        help='Name and BED file containing sites over which to restrict ' +
        'modified base results. Multiple sets of valid sites may be ' +
        'provided. For example E. coli 6mA sites could be specified as: ' +
        '`--valid-sites "Dam Methylation" Dam_motif_sites.bed ' +
        '--valid-sites "EcoKI Methylation" EcoKI_motif_sites.bed`.')
    mod_grp.add_argument(
        '--strand-specific-sites', action='store_true',
        help='Sites in --ground-truth-data and/or --valid-sites are ' +
        'strand-specific. Default: Sites are not strand specific.')
    mod_grp.add_argument(
        '--allow-unbalance-classes', action='store_true',
        help='Allow unbalanced classes in modified base metric computation. ' +
        'Default: Balance size of modified and canonical classes for each ' +
        'comparison made.')

    out_grp = parser.add_argument_group('Output Arguments')
    out_grp.add_argument(
        '--results-labels', nargs='+',
        help='Name for each Megalodon results directory. Control ' +
        'directories will have the suffix " Control" appended to the names. ' +
        'Default: "Sample 1", "Sample 2", ...')
    out_grp.add_argument(
        '--out-pdf', default='megalodon_validation.pdf',
        help='Output pdf filename. Default: %(default)s')
    out_grp.add_argument(
        '--out-filename',
        help='Output filename for text summary. Default: stdout')
    out_grp.add_argument(
        '--quiet', action='store_true',
        help='Suppress progress information.')

    return parser


def get_parser_validate_aggregated_modified_bases():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--modified-bed-methyl-files', nargs='+', required=True,
        help='Bed methyl files from modified sample(s).')
    parser.add_argument(
        '--ground-truth-csvs', nargs='+',
        help='Ground truth csvs with (chrm, strand, pos, is_mod) values. ' +
        'To collapse to forward strand coordinates, strand should be ".".')
    parser.add_argument(
        '--control-bed-methyl-files', nargs='+',
        help='Bed methyl files from control sample(s).')
    parser.add_argument(
        '--valid-positions', action='append',
        help='BED file containing positions to be considered. Multiple ' +
        'files may be provided')
    parser.add_argument(
        '--coverage-threshold', type=int, default=1,
        help='Only include sites with sufficient coverage. ' +
        'Default: 1 (= All sites)')
    parser.add_argument(
        '--strand-offset', type=int,
        help='Offset to combine stranded results. Positive value indicates ' +
        'reverse strand sites have higher position values. Default treat ' +
        'strands independently.')
    parser.add_argument(
        '--allow-unbalance-classes', action='store_true',
        help='Allow unbalanced classes in modified base metric computation. ' +
        'Default: Balance size of modified and canonical classes for each ' +
        'comparison made.')
    parser.add_argument(
        '--out-pdf', default='megalodon_agg_validation.pdf',
        help='Output pdf filename. Default: %(default)s')
    parser.add_argument(
        '--out-filename',
        help='Output filename for text summary. Default: stdout')

    return parser


def get_parser_validate_compare_modified_bases():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--sample1-bed-methyl-files', nargs='+', required=True,
        help='Bed methyl files from first set of sample(s).')
    parser.add_argument(
        '--sample2-bed-methyl-files', nargs='+', required=True,
        help='Bed methyl files from second set of sample(s).')
    parser.add_argument(
        '--sample-names', nargs=2, default=['Sample 1', 'Sample 2'],
        help='Name for provided samples. Default: %(default)s')
    parser.add_argument(
        '--valid-positions', action='append',
        help='BED file containing positions to be considered. Multiple ' +
        'files may be provided')
    parser.add_argument(
        '--coverage-threshold', type=int, default=1,
        help='Only include sites with sufficient coverage. ' +
        'Default: 1 (= All sites)')
    parser.add_argument(
        '--heatmap-num-bins', type=int, default=31,
        help='Number of bins for heatmap plotting. ' +
        'Default: %(default)d')
    parser.add_argument(
        '--strand-offset', type=int,
        help='Offset to combine stranded results. Positive value indicates ' +
        'reverse strand sites have higher position values. Default treat ' +
        'strands independently.')
    parser.add_argument(
        '--out-pdf', default='megalodon_mod_comaparison.pdf',
        help='Output pdf filename. Default: %(default)s')
    parser.add_argument(
        '--out-filename',
        help='Output filename for text summary. Default: stdout')

    return parser


###################
# Variant Parsers #
###################

def get_parser_variants_atomize():
    parser = argparse.ArgumentParser(
        description='Atomize variants so this does not have to be ' +
        'completed during read processing')
    parser.add_argument('in_vcf', help='Proposed varitants (VCF)')
    parser.add_argument(
        'reference',
        help='Reference FASTA or minimap2 index file corresponding to VCF.')
    parser.add_argument(
        '--out-vcf', default='atomized_variants.megalodon.vcf',
        help='Output VCF file. Default: %(default)s')
    parser.add_argument(
        '--max-indel-size', type=int, default=50,
        help='Maximum difference in number of reference and alternate ' +
        'bases. Default: %(default)d')
    return parser


def get_parser_variants_resolve():
    parser = argparse.ArgumentParser(
        description='Consolidate variants including filtering out ' +
        'reference variants and calling overlapping variants.')
    parser.add_argument(
        'variants',
        help='Megalodon called variant file. Must contain GL sample field.')
    parser.add_argument(
        '--output-filename', default='megalodon.consolidated_variants.vcf',
        help='Output filename. Default: %(default)s')
    parser.add_argument(
        '--max-likelihood-ratio', type=float, default=1,
        help='Maximum likelihood ratio ([ref prob] / [max alt prob]) to ' +
        'include variant in output. Allows output of uncertain reference ' +
        'calls. Default: 1; Include only sites called as alternative.')
    parser.add_argument(
        '--min-depth', type=int,
        help='Minimum depth to include a variant. Default: No depth filter')
    parser.add_argument(
        '--trim-variants', action='store_true',
        help='Trim extra padding sequence included by megalodon (e.g. ' +
        'around repeat-region indels). Default: Output as found in input ' +
        'variants.')

    ssv_grp = parser.add_argument_group('Strand-specific Variant Arguments')
    ssv_grp.add_argument(
        '--reverse-strand-variants',
        help='Variants file produced only from reads mapping to the reverse ' +
        'strand. If provided, this assumes that the main variants file ' +
        'contains variants only supported by reads from the forward strand. ' +
        'This is used to identify systematic basecalling error variants. ' +
        'Errors made on both strands indicate potential putative variants ' +
        'and are thus excluded. Homopolymer variants occuring on both ' +
        'strands are included by default. Exclude these variants as well ' +
        'by setting --exclude-both-strand-homopolymers .')
    ssv_grp.add_argument(
        '--homopolymer-min-length', type=int, default=4,
        help='Minimum length to consider a variant as a homopolymer. ' +
        'Default: %(default)d')
    ssv_grp.add_argument(
        '--exclude-both-strand-homopolymers', action='store_true',
        help='By default homopolymer variants are included even if they ' +
        'occur on both strands. Set this flag to treat homopolymer variants ' +
        'as other variants.')

    return parser


def get_parser_variants_heterozygous_factor():
    parser = argparse.ArgumentParser(
        description="""
        Given ground truth variants ground_truth.vcf and per_read_snp_calls.db
        from completed validation run:
        Example command line het testing:

        snp_h_fact=0.85
        indel_h_fact=0.78
        mkdir -p het_factor.$snp_h_fact.$indel_h_fact
        cp per_read_snp_calls.db het_factor.$snp_h_fact.$indel_h_fact/
        megalodon_extras aggregate run \
            --output-directory het_factor.$snp_h_fact.$indel_h_fact/ \
            --outputs snps --heterozygous-factor $snp_h_fact $indel_h_fact \
            --processes 8 --write-vcf-log-prob --reference reference.fa
        megalodon_extras variants heterozygous_factor \
            ground_truth.vcf het_factor.$snp_h_fact.$indel_h_fact/variants.vcf
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        'ground_truth_variants',
        help='VCF file containing ground truth diploid variant calls.')
    parser.add_argument(
        'megalodon_variants', default='megalodon_results/variants.vcf',
        help='VCF file containing diploid variant calls from megalodon.')

    return parser


def get_parser_variants_index_database():
    parser = argparse.ArgumentParser(
        description='Create per-read variant calls database index. Can ' +
        'rescue results from unexpected program crashes.')

    parser.add_argument(
        '--megalodon-directory', default='megalodon_results',
        help='Megalodon output directory containing per-read variant ' +
        'database to be indexed. Default: %(default)s')
    parser.add_argument(
        '--output-suffix', default='var_index_database',
        help='Log file suffix. Default: %(default)s')

    return parser


# all megalodon_extras command groups
GRP_AGG = 'aggregate'
CMD_AGG_RUN = 'run'

GRP_CALIB = 'calibrate'
CMD_CALIB_MODS = 'modified_bases'
CMD_CALIB_VARS = 'variants'
CMD_CALIB_GEN_MODS = 'generate_modified_base_stats'
CMD_CALIB_GEN_VARS = 'generate_variant_stats'
CMD_CALIB_MERGE_MODS = 'merge_modified_bases'

GRP_MERGE = 'merge'
CMD_MERGE_MODS = 'modified_bases'
CMD_MERGE_AGG_MODS = 'aggregated_modified_bases'
CMD_MERGE_VARS = 'variants'

GRP_MODS = 'modified_bases'
CMD_MODS_ALPHABET = 'describe_alphabet'
CMD_MODS_EST_THRESH = 'estimate_threshold'
CMD_MODS_UPDATE_DB = 'update_database'
CMD_MODS_GT = 'create_ground_truth'
CMD_MODS_MOTIF = 'create_motif_bed'
CMD_MODS_PER_SITE = 'per_site_thresholds'
CMD_MODS_INDEX = 'index_database'
CMD_MODS_SPLIT = 'split_by_motif'

GRP_PHASE = 'phase_variants'
CMD_PHASE_FILT_WHATSHAP = 'whatshap_filter'
CMD_PHASE_GET_HAP_READS = 'extract_haplotype_reads'
CMD_PHASE_MERGE_HAP = 'merge_haploid_variants'

GRP_TXT = 'per_read_text'
CMD_TXT_MODS = 'modified_bases'
CMD_TXT_VARS = 'variants'

GRP_VAL = 'validate'
CMD_VAL_RES = 'results'
CMD_VAL_AGG_MODS = 'aggregated_modified_bases'
CMD_VAL_COMP_MODS = 'compare_modified_bases'

GRP_VARS = 'variants'
CMD_VAR_ATOM = 'atomize'
CMD_VAR_RESOLVE = 'resolve'
CMD_VAR_HET_FACTOR = 'heterozygous_factor'
CMD_VAR_INDEX = 'index_database'

PARSERS = {
    GRP_AGG: {
        CMD_AGG_RUN: get_parser_aggregate_run},
    GRP_CALIB: {
        CMD_CALIB_MODS: get_parser_calibrate_modified_bases,
        CMD_CALIB_VARS: get_parser_calibrate_variants,
        CMD_CALIB_GEN_MODS: get_parser_calibrate_generate_modified_bases_stats,
        CMD_CALIB_GEN_VARS: get_parser_calibrate_generate_variants_stats,
        CMD_CALIB_MERGE_MODS: get_parser_calibrate_merge_modified_bases},
    GRP_MERGE: {
        CMD_MERGE_MODS: get_parser_merge_modified_bases,
        CMD_MERGE_AGG_MODS: get_parser_merge_aggregated_modified_bases,
        CMD_MERGE_VARS: get_parser_merge_variants},
    GRP_MODS: {
        CMD_MODS_ALPHABET: get_parser_modified_bases_describe_alphabet,
        CMD_MODS_EST_THRESH: get_parser_modified_bases_estimate_threshold,
        CMD_MODS_UPDATE_DB: get_parser_modified_bases_update_database,
        CMD_MODS_GT: get_parser_modified_bases_create_ground_truth,
        CMD_MODS_MOTIF: get_parser_modified_bases_create_motif_bed,
        CMD_MODS_PER_SITE: get_parser_modified_bases_per_site_thresholds,
        CMD_MODS_INDEX: get_parser_modified_bases_index_database,
        CMD_MODS_SPLIT: get_parser_modified_bases_split_calls_by_motif},
    GRP_PHASE: {
        CMD_PHASE_FILT_WHATSHAP: get_parser_phase_variants_whatshap_filter,
        CMD_PHASE_GET_HAP_READS:
        get_parser_phase_variants_extract_haplotype_reads,
        CMD_PHASE_MERGE_HAP: get_parser_phase_variants_merge_haploid_variants},
    GRP_TXT: {
        CMD_TXT_MODS: get_parser_per_read_text_modified_bases,
        CMD_TXT_VARS: get_parser_per_read_text_variants},
    GRP_VAL: {
        CMD_VAL_RES: get_parser_validate_results,
        CMD_VAL_AGG_MODS: get_parser_validate_aggregated_modified_bases,
        CMD_VAL_COMP_MODS: get_parser_validate_compare_modified_bases},
    GRP_VARS: {
        CMD_VAR_ATOM: get_parser_variants_atomize,
        CMD_VAR_RESOLVE: get_parser_variants_resolve,
        CMD_VAR_HET_FACTOR: get_parser_variants_heterozygous_factor,
        CMD_VAR_INDEX: get_parser_variants_index_database}}


if __name__ == '__main__':
    sys.stderr.write(
        'This is a module. See commands with `megalodon_extras -h`')
    sys.exit(1)
