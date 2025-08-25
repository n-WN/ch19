"""Educational Coppersmith toolkit (teaching-only).

Modules:
- poly: integer polynomial utilities (dict-based)
- lll: Fraction-based LLL reduction for integer matrices
- univariate: Howgrave–Graham-style univariate small-root search
- bivar / bivariate: bivariate poly ops and small-root search with elimination
- elimination: Bareiss determinant, Sylvester matrix, interpolation resultant

Notes:
- This package is designed for clarity and reproducibility, not speed or hardening.
- All APIs intentionally use standard library types and explicit integer arithmetic.
"""

from . import bivar, bivariate, elimination, lll, poly, univariate

__all__ = [
    "bivar",
    "bivariate",
    "elimination",
    "lll",
    "poly",
    "univariate",
]

__version__ = "0.1.0"
