# -*- coding: utf-8 -*-

""" Logger. """

import logging
import os

loggers = {}

def get_logger(name: str,
               out_path: str,
               ) -> logging.Logger:
    """ Get logger. """

    # 2023-04-25: making loggers as a global variable to prevent making multiple loggers in the GUI
    global loggers
    if loggers.get(name):
        return loggers.get(name)

    file_formatter = logging.Formatter('%(asctime)s~%(levelname)s~%(message)s~module:%(module)s~function:%(module)s')
    console_formatter = logging.Formatter('%(levelname)s -- %(message)s')

    file_handler = logging.FileHandler(os.path.join(out_path, "logfile.log"))
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(file_formatter)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)
    console_handler.setFormatter(console_formatter)

    logger = logging.getLogger(name)
    logger.propagate = False
    logger.setLevel(logging.INFO)

    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger

# def remove_logger(logger: logging.Logger) -> None:
#     """ Remove logger. """
#
#     handlers = logger.handlers[:]
#     for handler in handlers:
#         handler.close()
#         logger.removeHandler(handler)


