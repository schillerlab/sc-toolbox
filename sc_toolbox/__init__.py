"""Top-level package for sc-toolbox."""

__author__ = "Lukas Heumos"
__email__ = "lukas.heumos@posteo.net"
__version__ = "0.12.1"

from pypi_latest import PypiLatest

sc_toolbox_pypi_latest = PypiLatest("sc-toolbox", __version__)
sc_toolbox_pypi_latest.check_latest()

from . import plot as pl
from . import preprocessing as pp
from . import tools as tl
