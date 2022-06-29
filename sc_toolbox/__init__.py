"""Top-level package for sc-toolbox."""

__author__ = "Lukas Heumos"
__email__ = "lukas.heumos@posteo.net"
__version__ = "0.12.0"

from pypi_latest import PypiLatest

pertpy_pypi_latest = PypiLatest("sc-toolbox", __version__)
pertpy_pypi_latest.check_latest()

from . import plot as pl
from . import preprocessing as pp
from . import tools as tl
