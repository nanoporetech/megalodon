import sys
import argparse
from importlib import import_module

from megalodon import _version
from megalodon_extras import _extras_parsers as ep


NESTED_COMMANDS = [
    (ep.GRP_AGG,
     'Aggregate per-read modified base and/or variant statistics',
     [(ep.CMD_AGG_RUN,
       'Run aggregation'), ]),

    (ep.GRP_CALIB,
     'Calibrate model statistics with ground truth modified base or variants',
     [(ep.CMD_CALIB_MODS,
       'Calibrate modified base statistics'),
      (ep.CMD_CALIB_VARS,
       'Calibrate sequence variant statistics'),
      (ep.CMD_CALIB_GEN_MODS,
       'Generate statistics for modified base calibration'),
      (ep.CMD_CALIB_GEN_VARS,
       'Generate statistics for sequence variant calibration'),
      (ep.CMD_CALIB_MERGE_MODS,
       'Merge modified base calibration files')]),

    (ep.GRP_MERGE,
     'Merge per-read databases or aggregated files',
     [(ep.CMD_MERGE_MODS,
       'Merge per-read modified base database'),
      (ep.CMD_MERGE_AGG_MODS,
       'Merge aggregated modified base bedmethyl files'),
      (ep.CMD_MERGE_VARS,
       'Merge per-read sequence variants database')]),

    (ep.GRP_MODS,
     'Miscellaneous modified base operations',
     [(ep.CMD_MODS_ALPHABET,
       'Print the alphabet for a choosen model'),
      (ep.CMD_MODS_EST_THRESH,
       'Estimate optimal global modified base threshold for sequence markup'),
      (ep.CMD_MODS_UPDATE_DB,
       'Update modified base database from older versions of Megalodon'),
      (ep.CMD_MODS_GT,
       'Create ground truth modified base file from bedmethyl files'),
      (ep.CMD_MODS_MOTIF,
       'Create BED file of motif sites'),
      (ep.CMD_MODS_PER_SITE,
       'Extract per-site modified base thresholds for signal mapping ' +
       'sequence markup'),
      (ep.CMD_MODS_INDEX,
       'Create per-read modified base database index'),
      (ep.CMD_MODS_SPLIT,
       'Split modified base database by motif')]),

    (ep.GRP_PHASE,
     'Phase variants',
     [(ep.CMD_PHASE_FILT_WHATSHAP,
       'Filter variants not compatible with whatshap'),
      (ep.CMD_PHASE_GET_HAP_READS,
       'Extract read ids from haplotypes determined by whatshap'),
      (ep.CMD_PHASE_MERGE_HAP,
       'Merge variants from haploid calls')]),

    (ep.GRP_TXT,
     'Output per-read text files',
     [(ep.CMD_TXT_MODS,
       'Output per-read modified base statistics text file'),
      (ep.CMD_TXT_VARS,
       'Output per-read sequence variant statistics text file')]),

    (ep.GRP_VAL,
     'Validate per-read mapping and modified base results',
     [(ep.CMD_VAL_RES,
       'Validate per-read mappings and modified bases (if available)'),
      (ep.CMD_VAL_AGG_MODS,
       'Validate aggregated modified bases results'),
      (ep.CMD_VAL_COMP_MODS,
       'Compare aggregated modified base results (bedMethyl files')]),

    (ep.GRP_VARS,
     'Miscellaneous sequence variant operations',
     [(ep.CMD_VAR_ATOM,
       'Atomize variants for faster processing'),
      (ep.CMD_VAR_RESOLVE,
       'Resolve potentially conflicting variants'),
      (ep.CMD_VAR_HET_FACTOR,
       'Estimate optimal heterozygous factors for diploid variant calling'),
      # TODO variant database API does not allow opening for writing once
      # database exists.
      (ep.CMD_VAR_INDEX,
       '***** Stub for future implementation *****')])
]


class SubcommandHelpFormatter(argparse.RawDescriptionHelpFormatter):
    def _format_action(self, action):
        parts = super(SubcommandHelpFormatter, self)._format_action(action)
        if action.nargs == argparse.PARSER:
            parts = "\n".join(parts.split("\n")[1:])
        return parts


def _main():
    """ The main routine.
    """
    desc = ('Megalodon extras command groups (additional help available ' +
            'within each command group):\n' + '\n'.join([
                '\t{0: <25}{1}'.format(grp_name, grp_help)
                for grp_name, grp_help, _ in NESTED_COMMANDS]))
    parser = argparse.ArgumentParser(
        prog='megalodon_extras',
        description='********** Megalodon Extras *********\n\n' +
        'Commands to perform operations related to main Megalodon command ' +
        'including aggregation, variant phasing, validation, and more.\n\n' +
        desc, formatter_class=SubcommandHelpFormatter)
    parser.add_argument(
        '-v', '--version', action='version',
        version='Megalodon version: {}'.format(
            _version.MEGALODON_VERSION),
        help='Show Megalodon version and exit.')

    # add megalodon_extras command groups
    # service_command is a groupings of other commands
    # action_command is an executable command with detailed argument help
    service_subparsers = parser.add_subparsers(dest="service_command")
    for grp_name, grp_help, grp_sub_cmds in NESTED_COMMANDS:
        grp_desc = '\n'.join([
            '\t{0: <30}{1}'.format(cmd_name, cmd_help)
            for cmd_name, cmd_help in grp_sub_cmds])
        grp_parser = service_subparsers.add_parser(
            grp_name, formatter_class=SubcommandHelpFormatter,
            description=grp_desc)
        grp_subparser = grp_parser.add_subparsers(
            title=grp_name, dest="action_command")
        for cmd_name, cmd_help in grp_sub_cmds:
            # add each action parser to this service parser group
            grp_subparser.add_parser(
                cmd_name, add_help=False,
                parents=[ep.PARSERS[grp_name][cmd_name](), ])

    args = parser.parse_args()
    # if no service parser was provided print help and return
    if args.service_command is None:
        sys.stderr.write(
            '********** Please provide a megalodon_extras command group ' +
            'for further help. **********\n')
        parser.print_help()
        sys.exit(2)

    # if no action parser is provided print that command groups help
    if args.action_command is None:
        sys.stderr.write(
            '********** Please provide a command for further help. ' +
            '**********\n')
        parser.parse_args([args.service_command, '-h'])

    module = import_module('.{}_{}'.format(
        args.service_command, args.action_command), 'megalodon_extras')
    module._main(args)


if __name__ == '__main__':
    _main()
