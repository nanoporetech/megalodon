import os
import logging

from megalodon import megalodon_helper as mh

class CustomFormatter(logging.Formatter):
    err_fmt  = "*" * 100 + "\n\tERROR: %(msg)s\n" + "*" * 100
    warn_fmt  = "*" * 20 + " WARNING: %(msg)s " + "*" * 20
    info_fmt = "[%(asctime)s] %(message)s"
    dbg_fmt  = "DBG: %(module)s: %(lineno)d: %(msg)s"

    def __init__(self, fmt='[%(asctime)s] %(levelname)-8s: %(message)s'):
        super().__init__(fmt=fmt, datefmt='%H:%M:%S', style='%')
        return

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

def init_logger(out_dir, out_suffix=None):
    log_fn = os.path.join(out_dir, mh.LOG_FILENAME)
    if out_suffix is not None:
        base_fn, fn_ext = os.path.splitext(log_fn)
        log_fn = base_fn + '.' + out_suffix + fn_ext
    log_file = logging.FileHandler(log_fn, 'w')
    log_file.setLevel(logging.DEBUG)
    log_file.setFormatter(CustomFormatter())
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(CustomFormatter())

    root_logger = logging.getLogger('')
    root_logger.setLevel(logging.DEBUG)
    root_logger.addHandler(log_file)
    root_logger.addHandler(console)

    return

def get_logger(module_name=''):
    return logging.getLogger(module_name)


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
