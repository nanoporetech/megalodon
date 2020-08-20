import sys
from megalodon import logging, mods, megalodon_helper as mh
from ._extras_parsers import get_parser_modified_bases_index_database


LOGGER = logging.get_logger()


def _main(args):
    logging.init_logger(
        args.megalodon_directory, out_suffix=args.output_suffix)
    LOGGER.debug('Command: """' + ' '.join(sys.argv) + '"""')

    mods_db_fn = mh.get_megalodon_fn(args.megalodon_directory, mh.PR_MOD_NAME)
    mods_db = mods.ModsDb(mods_db_fn, read_only=False)
    try:
        mods_db.check_data_covering_index_exists()
        LOGGER.info('Modified bases database index already exists')
    except mh.MegaError:
        LOGGER.info('Creating modified bases database index')
        mods_db.create_data_covering_index()
    LOGGER.debug('Closing database')
    mods_db.close()


if __name__ == '__main__':
    _main(get_parser_modified_bases_index_database().parse_args())
