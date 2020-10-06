import argparse

from megalodon import backends, logging, megalodon_helper as mh


LOGGER = logging.get_logger()


def open_pyguppy_backend(args):
    args.do_not_use_guppy_server = False
    try:
        mh.mkdir(args.output_directory, False)
    except mh.MegaError:
        LOGGER.warning(
            'Guppy logs output directory exists. Potentially overwriting ' +
            'guppy logs.')
    backend_params = backends.parse_backend_params(args)
    model_info = None
    try:
        model_info = backends.ModelInfo(backend_params, args.processes)
        # if spawning multiple workers run this inside newly spawned processes
        model_info.prep_model_worker()
        LOGGER.info(model_info.get_alphabet_str())
        LOGGER.info('Model structure:\n\tStride: {}\n\tState size: {}'.format(
            model_info.stride, model_info.output_size))
        model_info.client.disconnect()
    finally:
        # ensure guppy server is closed in finally block
        if model_info is not None:
            model_info.close()


def get_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--log-directory', default='.',
        help='Directory to output megalodon log. Default: current ' +
        'working directory')

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
        '--output-directory', default='guppy_logs',
        help='Directory to output guppy logs. Default: %(default)s')
    pyg_grp.add_argument(
        '--devices', nargs='+',
        help='GPU devices for guppy basecalling backend.')
    pyg_grp.add_argument(
        '--processes', type=int, default=1,
        help='Number of parallel processes. Default: %(default)d')

    return parser


def main(args):
    logging.init_logger(args.log_directory)
    open_pyguppy_backend(args)


if __name__ == '__main__':
    main(get_parser().parse_args())
