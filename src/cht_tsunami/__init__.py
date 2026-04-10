"""cht_tsunami — tsunami source generation and initial condition utilities.

Exposes fault database access helpers and the :class:`Tsunami` class for
computing Okada (1985) sea-floor displacement fields.
"""

from .faults import (  # noqa: F401
    get_faults,
    get_okada_params_from_fault,
    parse_gem_tuple,
)
from .tsunami import Tsunami  # noqa: F401
