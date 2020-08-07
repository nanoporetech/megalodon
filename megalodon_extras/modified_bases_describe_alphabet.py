from megalodon import backends, logging, megalodon_helper as mh
from ._extras_parsers import get_parser_modified_bases_describe_alphabet


LOGGER = logging.get_logger()


def _main(args):
    logging.init_logger()
    # set args that are not relevant to alphabet
    args.devices = None

    # set guppy args
    args.guppy_server_port = None
    args.guppy_timeout = mh.DEFAULT_GUPPY_TIMEOUT
    args.output_directory = args.guppy_logs_output_directory

    # set taiyaki args
    args.chunk_size = 1000
    args.chunk_overlap = 100
    args.max_concurrent_chunks = 200
    try:
        mh.mkdir(args.output_directory, False)
    except mh.MegaError:
        LOGGER.warning(
            'Guppy logs output directory exists. Potentially overwriting ' +
            'guppy logs.')
    backend_params = backends.parse_backend_params(args)
    with backends.ModelInfo(backend_params, 1) as model_info:
        if model_info.is_cat_mod:
            LOGGER.info(
                'Using canonical alphabet {} and modified bases {}.'.format(
                    model_info.can_alphabet, '; '.join(
                        '{}={} (alt to {})'.format(
                            mod_b, mln, model_info.mod_base_to_can[mod_b])
                        for mod_b, mln in model_info.mod_long_names)))
        else:
            LOGGER.info(
                'Model contains canonical alphabet {}.'.format(
                    model_info.can_alphabet))


if __name__ == '__main__':
    _main(get_parser_modified_bases_describe_alphabet().parse_args())
