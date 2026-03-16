import sys
import logging
from .term_colors import TermColors as Col


def setup_log(name, debug=False):
    if "__main__" == name:
        name = "main"
    py_logger = logging.getLogger(name)
    log_level = logging.INFO if not debug else logging.DEBUG
    py_logger.setLevel(log_level)
    # py_formatter = logging.Formatter('tumorlens::%(asctime)s::%(levelname)s::%(name)s>  %(message)s')
    if not py_logger.handlers:
        py_handler = logging.StreamHandler(sys.stderr)
        py_handler.setFormatter(CustomFormatter())
        py_logger.addHandler(py_handler)
    return py_logger


class CustomFormatter(logging.Formatter):
    CRT = f'{Col.Fg.red}{Col.bold}'
    ERR = f'{Col.Fg.red}'
    DEB = f'{Col.Fg.green}'
    WRN = f'{Col.Fg.orange}'
    END = f'{Col.end}'
    format = 'tumorlens::%(asctime)s::%(levelname)s::%(name)s>  %(message)s'

    FORMATS = {
        logging.INFO: format,
        logging.DEBUG: DEB + format + END,
        logging.WARNING: WRN + format + END,
        logging.ERROR: ERR + format + END,
        logging.CRITICAL: CRT + format + END
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)
