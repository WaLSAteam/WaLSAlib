"""
WaLSAlib: WaLSA community library (auxiliary routines complementing WaLSAtools).

Convenience:
    from WaLSAlib import linefit
Optionally (if WaLSAtools is installed):
    from WaLSAlib import WaLSAtools
"""

# Export linefit as a top-level callable
from .linefit import linefit

__all__ = ["linefit"]

# Optional re-export of WaLSAtools (only if installed)
try:
    import WaLSAtools as WaLSAtools  # case-sensitive import
    __all__.append("WaLSAtools")
except Exception:
    # WaLSAtools not installed (or not importable in this environment)
    pass