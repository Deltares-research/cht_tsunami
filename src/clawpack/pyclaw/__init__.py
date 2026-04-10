"""Main pyclaw package"""

import logging
import logging.config
import os

_init = os.path.abspath(__file__)
_root = os.path.dirname(os.path.dirname(os.path.dirname(_init)))
if os.path.isdir(_root):
    __path__.append(_root)

# Default logging configuration file
_DEFAULT_LOG_CONFIG_PATH = os.path.join(os.path.dirname(__file__), "log.config")
del os, _init, _root

# Setup loggers
logging.config.fileConfig(_DEFAULT_LOG_CONFIG_PATH, disable_existing_loggers=False)

__all__ = []

# Module imports
__all__.extend(
    ["Controller", "Dimension", "Patch", "Domain", "Solution", "State", "CFL", "plot"]
)
from clawpack.pyclaw import plot

from .cfl import CFL
from .controller import Controller
from .geometry import Dimension, Domain, Patch
from .solution import Solution
from .state import State

__all__.extend(
    [
        "ClawSolver1D",
        "ClawSolver2D",
        "ClawSolver3D",
        "SharpClawSolver1D",
        "SharpClawSolver2D",
        "SharpClawSolver3D",
    ]
)
# Sub-packages
from . import limiters
from .classic.solver import ClawSolver1D, ClawSolver2D, ClawSolver3D
from .limiters import *
from .sharpclaw.solver import (
    SharpClawSolver1D,
    SharpClawSolver2D,
    SharpClawSolver3D,
)

__all__.extend(limiters.__all__)

__all__.append("BC")
from .solver import BC

__all__.extend("IOTest")
from .tests.test_io import IOTest
