import sys
import logging

LOGGING_FMT = "%(name)-20s %(levelname)-7s @ %(asctime)s: %(message)s"
LOGGING_DATE_FMT = "%m/%d/%y %H:%M:%S"

def get_logger(name, level=None):
    log = logging.getLogger(name)
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter(fmt=LOGGING_FMT, datefmt=LOGGING_DATE_FMT))
    log.addHandler(handler)
    if level:
        log.setLevel(level)
    else:
        log.setLevel(logging.DEBUG)
    return log
