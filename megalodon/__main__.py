import sys
import argparse

from megalodon import megalodon_helper as mh
from megalodon._version import MEGALODON_VERSION


class SelectiveRawFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        # special splitlines command for options output for better readability
        if text.startswith('O|'):
            return text[2:].splitlines()
        # else use standard RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)


def get_parser():
    if '--list-supported-guppy-configs' in sys.argv:
        sys.stderr.write(mh.get_supported_configs_message())
        sys.exit()

    # hide more complex arguments for standard help output
    show_hidden_args = '--help-long' in sys.argv

    def hidden_help(help_msg):
        if not show_hidden_args:
            return argparse.SUPPRESS
        return help_msg

    parser = argparse.ArgumentParser(formatter_class=SelectiveRawFormatter)
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
        '--do-not-use-guppy-server', action='store_true',
        help=hidden_help('Use alternative basecalling backend. Either ' +
                         'FAST5 (default; requires --post_out when running ' +
                         'guppy) or taiyaki (set `--taiyaki-model-filename` ' +
                         'to use taiyaki backend).'))
    pyg_grp.add_argument(
        '--guppy-params',
        help=hidden_help('Extra guppy server parameters. Main purpose for ' +
                         'optimal performance based on compute environment. ' +
                         'Quote parameters to avoid them being parsed by ' +
                         'megalodon.'))
    pyg_grp.add_argument(
        '--guppy-server-port', type=int,
        help=hidden_help('Guppy server port. Default: Guppy auto'))
    pyg_grp.add_argument(
        '--reads-per-guppy-batch', type=int,
        default=mh.DEFAULT_GUPPY_BATCH_SIZE,
        help=hidden_help('Number of reads to send to guppy per batch within ' +
                         'each worker processes. Default: %(default)d'))
    pyg_grp.add_argument(
        '--guppy-timeout', type=float, default=mh.DEFAULT_GUPPY_TIMEOUT,
        help=hidden_help('Timeout (in seconds) to wait for guppy server to ' +
                         'call a batch of reads. --processes and ' +
                         '--reads-per-guppy-batch have a strong bearing on ' +
                         'the appropriate value. Default: %(default)f'))
    pyg_grp.add_argument(
        '--list-supported-guppy-configs', action='store_true',
        help=hidden_help('List guppy configs with sequence variant and ' +
                         '(if applicable) modified base support.'))

    out_grp = parser.add_argument_group('Output Arguments')
    out_grp.add_argument(
        '--live-processing', action='store_true',
        help='Process reads from a live sequencing run. The [fast5s_dir] ' +
        'must be the base MinKNOW output directory. Megalodon will continue ' +
        'searching for FAST5 files until the file starting with ' +
        '"final_summary" is found.')
    out_grp.add_argument(
        '--outputs', nargs='+',
        default=['basecalls', ], choices=tuple(mh.OUTPUT_DESCS.keys()),
        # note 'O|' triggers raw formatting for this option alone
        help='O|Desired output(s).\nOptions:\n' +
        '\n'.join(('\t{}: {}'.format(*out_desc)
                   for out_desc in mh.OUTPUT_DESCS.items())) +
        '\nDefault: %(default)s')
    out_grp.add_argument(
        '--output-directory',
        default='megalodon_results',
        help='Directory to store output results. Default: %(default)s')
    out_grp.add_argument(
        '--overwrite', action='store_true',
        help='Overwrite output directory if it exists.')

    out_grp.add_argument(
        '--basecalls-format', choices=mh.BC_OUT_FMTS,
        default=mh.BC_OUT_FMTS[0],
        help=hidden_help('Basecalls output format. Choices: {}'.format(
            ', '.join(mh.BC_OUT_FMTS))))
    out_grp.add_argument(
        '--num-reads', type=int,
        help=hidden_help('Number of reads to process. Default: All reads'))
    out_grp.add_argument(
        '--read-ids-filename',
        help=hidden_help('File containing read ids to process (one per ' +
                         'line). Default: All reads'))

    map_grp = parser.add_argument_group('Mapping Arguments')
    map_grp.add_argument(
        '--mappings-format', choices=mh.MAP_OUT_FMTS,
        default=mh.MAP_OUT_FMTS[0],
        help='Mappings output format. Choices: {}'.format(
            ', '.join(mh.MAP_OUT_FMTS)))
    map_grp.add_argument(
        '--reference',
        help='Reference FASTA or minimap2 index file used for mapping ' +
        'called reads.')

    map_grp.add_argument(
        '--cram-reference',
        help=hidden_help('FASTA reference file. If --reference is a ' +
                         'minimap2 index, the associated FASTA reference ' +
                         'needs to be provided for the CRAM mapping output ' +
                         'format.'))
    map_grp.add_argument(
        '--samtools-executable', default='samtools',
        help=hidden_help('Samtools executable or path. Default: %(default)s'))
    map_grp.add_argument(
        '--sort-mappings', action='store_true',
        help=hidden_help('Perform sorting and indexing of mapping output ' +
                         'files. This can take considerable time for larger ' +
                         'runs.'))

    var_grp = parser.add_argument_group('Sequence Variant Arguments')
    var_grp.add_argument(
        '--haploid', action='store_true',
        help='Compute variant aggregation for haploid genotypes. ' +
        'Default: diploid')
    var_grp.add_argument(
        '--variant-filename',
        help='Sequence variants to call for each read in VCF/BCF format ' +
        '(required for variant output).')

    var_grp.add_argument(
        '--context-min-alt-prob', type=float,
        default=mh.DEFAULT_CONTEXT_MIN_ALT_PROB,
        help=hidden_help('Minimum alternative alleles probability to ' +
                         'include variant in computation of nearby variants.' +
                         ' Default: %(default)f'))
    var_grp.add_argument(
        '--disable-variant-calibration', action='store_true',
        help=hidden_help('Use raw variant scores from the network. ' +
                         'Default: Calibrate score with ' +
                         '--variant-calibration-filename'))
    var_grp.add_argument(
        '--heterozygous-factors', type=float, nargs=2,
        default=[mh.DEFAULT_SNV_HET_FACTOR, mh.DEFAULT_INDEL_HET_FACTOR],
        help=hidden_help('Bayesian prior factor for snv and indel ' +
                         'heterozygous calls. Smaller values preference ' +
                         'heterozygous calls; Larger values perference ' +
                         'homozygous calls. Default: %(default)s'))
    var_grp.add_argument(
        '--max-indel-size', type=int, default=mh.DEFAULT_MAX_INDEL_SIZE,
        help=hidden_help('Maximum difference in number of reference and ' +
                         'alternate bases. Default: %(default)d'))
    var_grp.add_argument(
        '--variant-all-paths', action='store_true',
        help=hidden_help('Compute forwards algorithm all paths score. ' +
                         '(Default: Viterbi best-path score)'))
    var_grp.add_argument(
        '--variants-are-atomized', action='store_true',
        help=hidden_help('Input variants have been atomized (with ' +
                         '`megalodon_extras variants atomize` command). ' +
                         'This saves compute time, but has unpredictable ' +
                         'behavior if variants are not atomized.'))
    var_grp.add_argument(
        '--variant-calibration-filename',
        help=hidden_help('File containing emperical calibration for ' +
                         'variant scores. See `megalodon_extras calibrate ' +
                         'variants` command. Default: Load default ' +
                         'calibration for specified guppy config.'))
    var_grp.add_argument(
        '--variant-context-bases', type=int, nargs=2,
        default=mh.DEFAULT_VAR_CONTEXT_BASES,
        help=hidden_help('Context bases for single base variant and indel ' +
                         'calling. Default: %(default)s'))
    var_grp.add_argument(
        '--variant-locations-on-disk', action='store_true',
        help=hidden_help('Force sequence variant locations to be stored ' +
                         'only within on disk database table. This option ' +
                         'will reduce the RAM memory requirement, but may ' +
                         'drastically slow processing. Default: Store ' +
                         'locations in memory and on disk.'))
    var_grp.add_argument(
        '--write-variants-text', action='store_true',
        help=hidden_help('Write per-read sequence variant calls out to a ' +
                         'text file. Default: Only ouput to database.'))
    var_grp.add_argument(
        '--write-vcf-log-probs', action='store_true',
        help=hidden_help('Write per-read alt log probabilities out in ' +
                         'non-standard VCF field.'))

    mod_grp = parser.add_argument_group('Modified Base Arguments')
    mod_grp.add_argument(
        '--mod-motif', action="append", nargs=3,
        metavar=('BASE', 'MOTIF', 'REL_POSITION'),
        help='Restrict modified base calls to specified motif(s). Argument ' +
        'takes 3 values representing 1) the single letter modified base(s), ' +
        '2) sequence motif and 3) relative modified base position. Multiple ' +
        '--mod-motif arguments may be provided to a single command. For ' +
        'example to restrict to CpG sites use "--mod-motif Z CG 0".')

    mod_grp.add_argument(
        '--disable-mod-calibration', action='store_true',
        help=hidden_help('Use raw modified base scores from the network. ' +
                         'Default: Calibrate scores as described in ' +
                         '--mod-calibration-filename'))
    mod_grp.add_argument(
        '--mod-aggregate-method', choices=list(mh.MOD_AGG_METHOD_NAMES),
        default=mh.MOD_BIN_THRESH_NAME,
        help=hidden_help('Modified base aggregation method. ' +
                         'Default: %(default)s'))
    mod_grp.add_argument(
        '--mod-all-paths', action='store_true',
        help=hidden_help('Compute forwards algorithm all paths score for ' +
                         'modified base calls. (Default: Viterbi ' +
                         'best-path score)'))
    out_grp.add_argument(
        '--mod-min-prob', type=float, default=mh.DEFAULT_MOD_MIN_PROB,
        help=hidden_help('Only include modified base probabilities greater ' +
                         'than this value in mod_basecalls and mod_mappings ' +
                         'outputs. Default: %(default)f'))
    mod_grp.add_argument(
        '--mod-binary-threshold', type=float,
        default=mh.DEFAULT_MOD_BINARY_THRESH,
        help=hidden_help('Threshold for modified base aggregation ' +
                         '(probability of modified/canonical base). ' +
                         'Only applicable for "--mod-aggregate-method ' +
                         'binary_threshold". Default: %(default)s'))
    mod_grp.add_argument(
        '--mod-calibration-filename',
        help=hidden_help('File containing emperical calibration for ' +
                         'modified base scores. See `megalodon_extras ' +
                         'calibrate modified_bases` command. Default: ' +
                         'Load default calibration for specified guppy ' +
                         'config.'))
    mod_grp.add_argument(
        '--mod-database-timeout', type=float,
        default=mh.DEFAULT_MOD_DATABASE_TIMEOUT,
        help=hidden_help('Timeout in seconds for modified base database ' +
                         'operations. Default: %(default)f'))
    mod_grp.add_argument(
        '--mod-context-bases', type=int, default=mh.DEFAULT_MOD_CONTEXT,
        help=hidden_help('Context bases for modified base calling. ' +
                         'Default: %(default)d'))
    mod_grp.add_argument(
        '--mod-map-emulate-bisulfite', action='store_true',
        help=hidden_help('For mod_mappings output, emulate bisulfite output ' +
                         'by converting called bases setting ' +
                         '"--mod-map-base-conv" argument.'))
    mod_grp.add_argument(
        '--mod-map-base-conv', action='append', nargs=2,
        metavar=('FROM_BASE', 'TO_BASE'),
        help=hidden_help('For mod_mappings output, convert called modified ' +
                         'bases. Only applicable when ' +
                         '--mod-map-emulate-bisulfite is set.For example, ' +
                         'to emulate bisulfite output use: ' +
                         '"--mod-map-base-conv C T --mod-map-base-conv m C"'))
    mod_grp.add_argument(
        '--mod-output-formats', nargs='+',
        default=[mh.MOD_BEDMETHYL_NAME, ],
        choices=tuple(mh.MOD_OUTPUT_FMTS.keys()),
        help=hidden_help('Modified base aggregated output format(s). ' +
                         'Default: %(default)s'))
    mod_grp.add_argument(
        '--mod-per-site-threshold',
        help=hidden_help('BED file containing per-site thresholds for ' +
                         'marking up modified base references. ' +
                         'See scripts/per_site_markup.py'))
    mod_grp.add_argument(
        '--write-mod-log-probs', action='store_true',
        help=hidden_help('Write per-read modified base log probabilities ' +
                         'out in non-standard modVCF field.'))
    mod_grp.add_argument(
        '--write-mods-text', action='store_true',
        help=hidden_help('Write per-read modified bases out to a text ' +
                         'file. Default: Only ouput to database.'))

    tai_grp = parser.add_argument_group('Taiyaki Backend Arguments')
    tai_grp.add_argument(
        '--chunk-size', type=int, default=1000,
        help=hidden_help('Chunk length for base calling. ' +
                         'Default: %(default)d'))
    tai_grp.add_argument(
        '--chunk-overlap', type=int, default=100,
        help=hidden_help('Overlap between chunks to be stitched together. ' +
                         'Default: %(default)d'))
    tai_grp.add_argument(
        '--max-concurrent-chunks', type=int, default=200,
        help=hidden_help('Only process N chunks concurrently per-read (to ' +
                         'avoid GPU memory errors). Default: %(default)d'))
    tai_grp.add_argument(
        '--taiyaki-model-filename',
        help=hidden_help('Taiyaki basecalling model checkpoint file.'))

    sigmap_grp = parser.add_argument_group(
        'Reference/Signal Mapping Output Arguments')
    sigmap_grp.add_argument(
        '--ref-include-mods', action='store_true',
        help=hidden_help('Include modified base calls in signal_mappings/' +
                         'per_read_refs output.'))
    sigmap_grp.add_argument(
        '--ref-include-variants', action='store_true',
        help=hidden_help('Include variant calls in per_read_refs output ' +
                         '(does not apply to signal_mappings output).'))
    sigmap_grp.add_argument(
        '--ref-length-range', type=int, nargs=2,
        metavar=('MIN_LENGTH', 'MAX_LENGTH'),
        help=hidden_help('Only include reads with specified read length ' +
                         'in signal_mappings/per_read_refs output.'))
    sigmap_grp.add_argument(
        '--ref-percent-identity-threshold', type=float,
        help=hidden_help('Only include reads with higher percent identity ' +
                         'in signal_mappings/per_read_refs output.'))
    sigmap_grp.add_argument(
        '--ref-percent-coverage-threshold', type=float,
        help=hidden_help('Only include reads with higher read alignment ' +
                         'coverage in signal_mappings/per_read_refs output.'))
    sigmap_grp.add_argument(
        '--ref-mods-all-motifs', nargs=4, action='append',
        metavar=('MOD', 'MOD_LONG_NAME', 'MOTIF', 'REL_POS'),
        help=hidden_help('Annotate all motifs as modified (e.g. bacterial ' +
                         'methylase). This will ignore modified base calls ' +
                         'made and --mod-motif. Multiple mods and motifs ' +
                         'can be provided. Arguments should be structured ' +
                         'as --ref-mods-all-motifs [single letter modified ' +
                         'base code] [modified base long name] ' +
                         '[canonical motif] [mod position relative to motif]'))
    sigmap_grp.add_argument(
        '--ref-mod-threshold', type=float, default=0.0,
        help=hidden_help('Threshold (log(can_prob/mod_prob)) used to ' +
                         'annotate a modified bases in signal_mappings/' +
                         'per_read_refs output. See `megalodon_extras ' +
                         'modified_bases estimate_threshold` command. ' +
                         'Default: %(default)f'))

    mp_grp = parser.add_argument_group('Compute Resource Arguments')
    mp_grp.add_argument(
        '--processes', type=int, default=1,
        help='Number of parallel processes. Default: %(default)d')
    mp_grp.add_argument(
        '--devices', nargs='+',
        help='GPU devices for guppy or taiyaki basecalling backends.')

    mp_grp.add_argument(
        '--num-read-enumeration-threads', type=int,
        default=mh.DEFAULT_READ_ENUM_TS,
        help=hidden_help('Number of parallel threads to use for read ' +
                         'enumeration. Increase if input queue remains ' +
                         'empty, generally due to single read format FAST5s ' +
                         'or slow disk. Default: %(default)d'))
    mp_grp.add_argument(
        '--num-extract-signal-processes', type=int,
        default=mh.DEFAULT_EXTRACT_SIG_PROC,
        help=hidden_help('Number of parallel processes to use for signal ' +
                         'extraction. Increasing this value can allow more ' +
                         'efficient raw data extraction. Note that ' +
                         '[--num-read-enumeration-threads] will be opened ' +
                         'in each extract signal process. ' +
                         'Default: %(default)d'))

    misc_grp = parser.add_argument_group('Miscellaneous Arguments')
    misc_grp.add_argument(
        '--help-long', help='Show all options.', action='help')
    misc_grp.add_argument(
        '--rna', action='store_true',
        help='RNA input data. Requires RNA model. Default: DNA input data')
    misc_grp.add_argument(
        '-v', '--version', action='version',
        version='Megalodon version: {}'.format(MEGALODON_VERSION),
        help='show megalodon version and exit.')

    misc_grp.add_argument(
        '--database-safety', type=int, default=0,
        help=hidden_help('Setting for database performance versus ' +
                         'corruption protection. Options: 0 (DB corruption ' +
                         'on application crash), 1 (DB corruption on system ' +
                         'crash), 2 (DB safe mode). Default: %(default)d'))
    misc_grp.add_argument(
        '--skip-database-index', action='store_true',
        help=hidden_help('When outputting per_read_mods and not aggregated ' +
                         'mods output, skip database indexing. ' +
                         '"megalodon_extras modified_bases index_database" ' +
                         'must be run before downstream processing. This ' +
                         'can be useful to minimize time on GPU compute ' +
                         'resources. Will apply to variants output in ' +
                         'the future.'))
    misc_grp.add_argument(
        '--edge-buffer', type=int, default=mh.DEFAULT_EDGE_BUFFER,
        help=hidden_help('Do not process sequence variant or modified base ' +
                         'calls near edge of read mapping. ' +
                         'Default: %(default)d'))
    misc_grp.add_argument(
        '--not-recursive', action='store_true',
        help=hidden_help('Only search for fast5 read files directly found ' +
                         'within the fast5 directory. Default: search ' +
                         'recursively'))
    misc_grp.add_argument(
        '--suppress-progress-bars', action='store_true',
        help=hidden_help('Suppress progress bars output.'))
    misc_grp.add_argument(
        '--suppress-queues-status', action='store_true',
        help=hidden_help('Suppress dynamic status of output queues. Helpful ' +
                         'for diagnosing I/O issues.'))
    misc_grp.add_argument(
        '--verbose-read-progress', type=int, default=3,
        help=hidden_help('Output verbose output on read progress. Outputs ' +
                         'N most common points where reads could not be ' +
                         'processed further. Default: %(default)d'))

    return parser


def _main():
    args = get_parser().parse_args()
    # only import megalodon if actually processing (not just printing help)
    from megalodon import megalodon
    megalodon._main(args)


if __name__ == '__main__':
    _main()
