# -*- coding: utf-8 -*-

"""Exception hierarchy for Riana.

These replace the prior mix of ``assert``, ``sys.exit``, and generic
``Exception`` raises so that callers (CLI, GUI, library users) can
distinguish error categories programmatically.
"""


class RianaError(Exception):
    """Base exception for all Riana errors."""


class DataError(RianaError):
    """Input data missing, malformed, or inconsistent (files, PSMs, mzML)."""


class IntegrationError(RianaError):
    """Failure during isotopomer extraction or integration."""


class ModelingError(RianaError):
    """Failure during kinetic model fitting."""


class ValidationError(RianaError):
    """Invalid parameter, configuration, or argument value."""
