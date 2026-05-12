# -*- coding: utf-8 -*-

""" Logger. """

import logging
import os

# Cache loggers by (logger_name, output_dir) so concurrent runs writing to
# different directories do not share a file handler. The prior version cached
# by logger name alone, which silently sent later runs' logs to the first
# run's output directory.
_loggers: dict[tuple[str, str], logging.Logger] = {}


def get_logger(name: str,
               out_path: str,
               ) -> logging.Logger:
    """Return a logger that writes to ``out_path/logfile.log`` and to stderr.

    Cached per ``(name, out_path)`` so repeated calls within a process reuse
    handlers, but distinct output directories get distinct handlers.
    """

    key = (name, os.path.abspath(out_path))
    cached = _loggers.get(key)
    if cached is not None:
        return cached

    file_formatter = logging.Formatter(
        '%(asctime)s~%(levelname)s~%(message)s~module:%(module)s~function:%(funcName)s'
    )
    console_formatter = logging.Formatter('%(levelname)s -- %(message)s')

    file_handler = logging.FileHandler(os.path.join(out_path, "logfile.log"))
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(file_formatter)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)
    console_handler.setFormatter(console_formatter)

    logger = logging.getLogger(f'{name}::{key[1]}')
    logger.propagate = False
    logger.setLevel(logging.INFO)

    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    _loggers[key] = logger
    return logger
