import sys
from megalodon import logging, variants, megalodon_helper as mh
from ._extras_parsers import get_parser_variants_index_database


LOGGER = logging.get_logger()


def _main(args):
    raise NotImplementedError(
        'Variant index creation not currently implemented.')

    logging.init_logger(
        args.megalodon_directory, out_suffix=args.output_suffix)
    LOGGER.debug('Command: """' + ' '.join(sys.argv) + '"""')

    vars_db_fn = mh.get_megalodon_fn(args.megalodon_directory, mh.PR_VAR_NAME)
    vars_db = variants.VarsDb(vars_db_fn, read_only=False)
    try:
        vars_db.check_data_covering_index_exists()
        LOGGER.info('Variants database index already exists')
    except mh.MegaError:
        LOGGER.info('Creating variants database index')
        vars_db.create_data_covering_index()
    LOGGER.debug('Closing database')
    vars_db.close()


if __name__ == '__main__':
    _main(get_parser_variants_index_database().parse_args())
