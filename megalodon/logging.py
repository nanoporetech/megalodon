import os
import sys
import logging


LOG_FN = 'log.txt'


class CustomFormatter(logging.Formatter):
    err_fmt = "*" * 100 + "\n\tERROR: %(msg)s\n" + "*" * 100
    warn_fmt = "*" * 20 + " WARNING: %(msg)s " + "*" * 20
    info_fmt = "[%(asctime)s] %(message)s"
    dbg_fmt = ("DBG %(asctime)s : %(msg)s --- %(processName)s-" +
               "%(threadName)s %(module)s.py:%(lineno)d")

    def __init__(self, fmt='[%(asctime)s] %(levelname)-8s: %(message)s'):
        super().__init__(fmt=fmt, datefmt='%H:%M:%S', style='%')

    def format(self, record):
        format_orig = self._fmt

        # Replace the original format with one customized by logging level
        if record.levelno == logging.DEBUG:
            self._style._fmt = self.dbg_fmt
        elif record.levelno == logging.INFO:
            self._style._fmt = self.info_fmt
        elif record.levelno == logging.WARNING:
            self._style._fmt = self.warn_fmt
        elif record.levelno == logging.ERROR:
            self._style._fmt = self.err_fmt
        result = logging.Formatter.format(self, record)

        self._fmt = format_orig

        return result


def init_logger(out_dir=None, out_suffix=None, log_fn=None, quiet=False):
    """ Prepare logging output. Output file will be opened if out_dir or log_fn
    are specified. out_suffix will be added to the standard log.txt filename in
    out_dir (does not apply when log_fn is specified).

    File will include debug and above messages while stderr will include info
    and above. If quiet=True, stderr will include warning and above only.
    """
    log_fp = None
    if out_dir is not None:
        log_fn = os.path.join(out_dir, LOG_FN)
        if out_suffix is not None:
            base_fn, fn_ext = os.path.splitext(log_fn)
            log_fn = base_fn + '.' + out_suffix + fn_ext
    if log_fn is not None:
        log_fp = logging.FileHandler(log_fn, 'w')
        log_fp.setLevel(logging.DEBUG)
        log_fp.setFormatter(CustomFormatter())

    console = logging.StreamHandler()
    if quiet:
        console.setLevel(logging.WARNING)
    else:
        console.setLevel(logging.INFO)
    console.setFormatter(CustomFormatter())

    root_logger = logging.getLogger('')
    root_logger.setLevel(logging.DEBUG)
    if log_fp is not None:
        root_logger.addHandler(log_fp)
    root_logger.addHandler(console)


def get_logger(module_name=''):
    return logging.getLogger(module_name)


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
